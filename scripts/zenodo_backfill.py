#!/usr/bin/env python3
"""Back-fill Zenodo DOIs for GitHub releases published before the Zenodo
GitHub integration was enabled on the repository.

Zenodo only ingests releases created *after* a repository is toggled on in its
GitHub settings. This script replays the "release published" webhook payload to
Zenodo's GitHub events receiver for each pre-existing release, oldest to newest,
so Zenodo creates a record and mints a version DOI for each. Adapted from
imcf/zenodo-doi-for-existing-release, with the webhook token taken from the
environment (never committed) and a dry-run mode.

Prerequisites
  1. The repository is already enabled in Zenodo (Profile -> GitHub -> toggle
     ON). Enabling installs a webhook on the GitHub repo whose Payload URL ends
     in ?access_token=<TOKEN>.
  2. Export that token:  export ZENODO_WEBHOOK_TOKEN=<TOKEN>
     Find it at https://github.com/<owner>/<repo>/settings/hooks -> Edit the
     Zenodo hook -> copy everything after access_token= in the Payload URL.

Usage
  python scripts/zenodo_backfill.py --repo ALS-RSOXS/refloxide --dry-run
  python scripts/zenodo_backfill.py --repo ALS-RSOXS/refloxide

Caveats
  - Each non-dry-run POST publishes a Zenodo record and mints a permanent DOI.
    DOIs cannot be deleted. Run --dry-run first and confirm the release list.
  - Metadata for an old release is taken from the source archive of that tag.
    Tags without a .zenodo.json or CITATION.cff get Zenodo's default metadata
    (author = GitHub handle, etc.); edit each record's metadata in the Zenodo
    web UI afterward, or add the citation file before cutting future releases.
  - Targets the legacy receiver /api/hooks/receivers/github/events/. If Zenodo
    has retired it on the current platform the POST will not return 202/success;
    inspect the printed status and stop rather than retrying blindly.
"""

import argparse
import os
import sys
import time

import dotenv

dotenv.load_dotenv()

import requests

GITHUB_API = "https://api.github.com"
ZENODO_RECEIVER = "https://zenodo.org/api/hooks/receivers/github/events/"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo", required=True, help="owner/repo, e.g. ALS-RSOXS/refloxide"
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="list releases; do not POST to Zenodo"
    )
    parser.add_argument(
        "--sleep", type=float, default=3.0, help="seconds between submissions"
    )
    args = parser.parse_args()

    token = os.environ.get("ZENODO_WEBHOOK_TOKEN")
    if not args.dry_run and not token:
        print(
            "ZENODO_WEBHOOK_TOKEN is not set. Export it or use --dry-run.",
            file=sys.stderr,
        )
        return 2

    gh_headers = {"Accept": "application/vnd.github+json"}
    gh_pat = os.environ.get("GITHUB_TOKEN")
    if gh_pat:
        gh_headers["Authorization"] = f"Bearer {gh_pat}"

    repo_resp = requests.get(
        f"{GITHUB_API}/repos/{args.repo}", headers=gh_headers, timeout=30
    )
    repo_resp.raise_for_status()
    repo_obj = repo_resp.json()

    rel_resp = requests.get(
        f"{GITHUB_API}/repos/{args.repo}/releases",
        headers=gh_headers,
        params={"per_page": 100},
        timeout=30,
    )
    rel_resp.raise_for_status()
    releases = rel_resp.json()
    if not isinstance(releases, list) or not releases:
        print(f"No releases found for {args.repo}.")
        return 1

    ordered = list(
        reversed(releases)
    )  # GitHub returns newest first; publish oldest first
    print(f"{len(ordered)} release(s) for {args.repo}, oldest to newest:\n")

    for release in ordered:
        tag = release.get("tag_name")
        name = release.get("name") or tag
        published = release.get("published_at")
        print(f"  {name}  (tag {tag}, published {published})")
        if args.dry_run:
            continue

        payload = {"action": "published", "release": release, "repository": repo_obj}
        resp = requests.post(
            ZENODO_RECEIVER, params={"access_token": token}, json=payload, timeout=60
        )
        ok = resp.status_code in (200, 201, 202)
        print(f"    -> Zenodo HTTP {resp.status_code} {'OK' if ok else 'CHECK'}")
        if not ok:
            print(f"    response body: {resp.text[:500]}")
            print(
                "    Stopping. The legacy receiver may be unavailable or the token is wrong."
            )
            return 1
        time.sleep(args.sleep)

    if args.dry_run:
        print("\nDry run only. Re-run without --dry-run to mint DOIs (irreversible).")
    else:
        print(
            "\nDone. Check https://zenodo.org/ -> your uploads for the new records and DOIs."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
