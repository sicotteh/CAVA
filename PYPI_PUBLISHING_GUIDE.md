# Quick Guide: Publishing CAVA to PyPI

## Pre-Publication Checklist

- [ ] Git tag created: `git tag v2.0.14`
- [ ] CHANGELOG updated with release notes
- [ ] All tests pass: `pytest cava/test/`
- [ ] README is up-to-date
- [ ] No uncommitted changes: `git status`

## Step 1: Create PyPI Account

1. Go to https://pypi.org/account/register/
2. Create account and enable 2FA
3. Go to https://pypi.org/manage/account/
4. Create API token (scoped to project)
5. Save token securely

## Step 2: Test on TestPyPI (Recommended)

```bash
# Install build tools
pip install --upgrade build twine

# Navigate to project directory
cd /Users/m037385/Documents/AAAprojects/CAVA/CAVA

# Clean old builds
rm -rf dist/ build/ *.egg-info

# Build distribution
python -m build

# Validate
twine check dist/*

# Upload to TestPyPI
twine upload --repository testpypi dist/*
  # When prompted for username, use: __token__
  # When prompted for password, paste your TestPyPI token

# Test installation from TestPyPI
pip install --index-url https://test.pypi.org/simple/ CAVA
```

## Step 3: Publish to PyPI (Production)

Once TestPyPI validation is complete:

```bash
# Clean old builds
rm -rf dist/ build/ *.egg-info

# Build distribution
python -m build

# Validate again
twine check dist/*

# Upload to PyPI
twine upload dist/*
  # When prompted for username, use: __token__
  # When prompted for password, paste your PyPI token
```

## Step 4: Verify Publication

```bash
# Check package page
open https://pypi.org/project/CAVA/

# Install from PyPI
pip install CAVA

# Verify installation
python -c "import cava; print(f'CAVA {cava.__version__}')"
```

## Updating Version for Next Release

1. Update `cava/VERSION` file:
   ```
   2.0.15
   ```

2. Update CHANGELOG.md

3. Commit and tag:
   ```bash
   git add cava/VERSION CHANGELOG.md
   git commit -m "Bump version to 2.0.15"
   git tag v2.0.15
   git push origin main --tags
   ```

4. Build and publish (repeat Step 2/3 above)

## If You Need to Update a Release

Note: PyPI doesn't allow re-uploading the same version. You must:

1. Increment the patch version
2. Update `cava/VERSION`
3. Rebuild: `python -m build`
4. Re-upload: `twine upload dist/*`

## Troubleshooting

### "File already exists"
This happens if you try to upload the same version twice. Increment the version number.

### "Invalid distribution"
Run `twine check dist/*` to see validation errors.

### "401 Unauthorized"
Make sure you're using `__token__` as username and the correct token value as password.

### Long build time
Normal for packages with Cython extensions. Be patient! (~15-30 min)

## Environment Variables (Optional Setup)

For automated uploads, you can set credentials:

```bash
# Create ~/.pypirc file (be careful with permissions!)
cat > ~/.pypirc << EOF
[distutils]
index-servers =
    pypi
    testpypi

[pypi]
repository = https://upload.pypi.org/legacy/
username = __token__
password = pypi-AgEIcHlwaS5vcmc...

[testpypi]
repository = https://test.pypi.org/legacy/
username = __token__
password = pypi-AgEIcHlwaS5vcmc...
EOF

chmod 600 ~/.pypirc
```

Then upload with just:
```bash
twine upload dist/*
```

## Additional Resources

- PyPI: https://pypi.org/
- Test PyPI: https://test.pypi.org/
- Twine Docs: https://twine.readthedocs.io/
- Packaging Guide: https://packaging.python.org/

---

**Important:** Always test on TestPyPI first before uploading to production PyPI!

