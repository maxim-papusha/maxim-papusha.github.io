# Deployment Troubleshooting Guide

## Issue: Website Pages Not Deployed

### Root Cause

The website is not being deployed because **GitHub Pages is not enabled** in the repository settings.

### Evidence

The GitHub Actions workflow runs successfully through the build phase but fails during deployment with the following error:

```
Error: Failed to create deployment (status: 404) with build version 96df386003b5691782705ee2bfdb6fecbabda160. 
Ensure GitHub Pages has been enabled: https://github.com/maxim-papusha/maxim-papusha.github.io/settings/pages
```

**Workflow Run ID:** 21963275820 (failed on main branch on 2026-02-12)

### Solution

To fix this issue, GitHub Pages must be enabled and configured to use GitHub Actions as the deployment source:

1. **Go to Repository Settings**
   - Navigate to: https://github.com/maxim-papusha/maxim-papusha.github.io/settings/pages

2. **Configure GitHub Pages**
   - Under "Build and deployment"
   - Set **Source** to: `GitHub Actions`

3. **Trigger a New Deployment**
   - Push a commit to the `main` branch, or
   - Manually trigger the workflow from the Actions tab

### Verification

After enabling GitHub Pages, verify the deployment:

1. Check the GitHub Actions workflow runs at:
   https://github.com/maxim-papusha/maxim-papusha.github.io/actions

2. The deployment job should complete successfully

3. The website should be accessible at:
   https://maxim-papusha.github.io/

### Current Build Status

The MkDocs build process is working correctly:
- ✅ Dependencies install successfully
- ✅ `mkdocs build --strict` completes without errors
- ✅ Static site is generated in the `site/` directory
- ✅ Build job in GitHub Actions passes
- ❌ Deployment job fails due to GitHub Pages not being enabled

### Workflow Configuration

The deployment workflow (`.github/workflows/deploy.yml`) is correctly configured with:
- Build job that creates the site artifact
- Deploy job that deploys to GitHub Pages (requires Pages to be enabled)
- Proper permissions: `pages: write`, `id-token: write`
- Correct deployment action: `actions/deploy-pages@v4`

### Next Steps

**For Repository Maintainers:**
1. Enable GitHub Pages with "GitHub Actions" as the source
2. Verify deployment succeeds
3. Confirm website is accessible at the GitHub Pages URL

**Note:** This is a configuration issue, not a code issue. No code changes are required.
