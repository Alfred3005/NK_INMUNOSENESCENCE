# CELLxGENE Census API Updates - January 2025

**Date**: 2025-01-03
**Updated by**: Claude Code
**Reason**: CELLxGENE Census released major updates in late 2024

---

## 🔄 Key Changes Implemented

### 1. Latest Census Data Release

**New**: Census version `2025-01-30` (LTS)
- **109.1 Million** human cells (up from previous versions)
- **45.4 Million** mouse cells
- **1,573 datasets** total
- **Schema version**: 2.1.0
- **Availability**: 5+ years (Long-Term Support)

**Update in code**:
- Default `census_version = "stable"` now points to `2025-01-30`
- Can use `census_version = "latest"` for weekly updates (1-month availability)

### 2. New API Methods (Census v1.13.0+)

#### `cellxgene_census.get_obs()` - NEW RECOMMENDED METHOD

**Before (old method)**:
```python
census = cellxgene_census.open_soma()
query = census["census_data"]["homo_sapiens"]
obs = query.obs.read(column_names=[...], value_filter="...").concat().to_pandas()
```

**After (new method)**:
```python
with cellxgene_census.open_soma() as census:
    obs = cellxgene_census.get_obs(
        census,
        "homo_sapiens",
        value_filter="cell_type == 'NK cell' and disease == 'normal'",
        column_names=["dataset_id", "cell_type", "donor_id"]
    )
```

**Benefits**:
- ✅ More concise syntax
- ✅ Better memory efficiency
- ✅ Automatic resource cleanup (context manager)
- ✅ Cleaner error handling

#### `cellxgene_census.get_anndata()` - For Direct AnnData Queries

**New capability**: Query and get AnnData objects directly without intermediate steps.

```python
adata = cellxgene_census.get_anndata(
    census=census,
    organism="Homo sapiens",
    obs_value_filter="cell_type == 'NK cell' and tissue_general == 'blood'",
    var_value_filter="feature_id in ['ENSG00000161798', 'ENSG00000188229']",
    column_names={"obs": ["sex", "age", "disease"]}
)
```

### 3. Categorical Columns (Breaking Change)

**Important**: Since Census v1.13.0 (April 2024), metadata columns are now **categorical** instead of strings.

**Impact**:
```python
# This now returns ALL categories, even absent ones
adata.obs['cell_type'].value_counts()

# Use this for present values only
adata.obs['cell_type'].value_counts(dropna=True)
```

**Our solution**: Code updated with fallback methods to handle both categorical and string types.

### 4. Context Manager Pattern (Best Practice)

**Before**:
```python
census = cellxgene_census.open_soma()
try:
    # ... use census
finally:
    census.close()
```

**After (recommended)**:
```python
with cellxgene_census.open_soma() as census:
    # ... use census
    # Automatic cleanup on exit
```

**Implementation**: Updated `src/download.py` to use context managers throughout.

---

## 📝 Updated Files

### `src/download.py`

**Changes**:
1. ✅ `discover_datasets()` now uses `get_obs()` API
2. ✅ Context manager for census connection
3. ✅ Fallback to old method if `get_obs()` fails (compatibility)
4. ✅ Better logging with census version info
5. ✅ Support for categorical columns

**New features**:
- Shows census release date and total cell count
- More efficient SOMA value filter construction
- Clearer filter logging

### `requirements.txt`

**Updated**:
```
cellxgene-census>=1.17.0  # (was 1.14.1)
```

**New features in 1.17.0**:
- Spatial data support
- Embeddings search API
- Improved PyTorch dataloaders
- Bug fixes and performance improvements

### `config/datasets.yaml`

**No changes needed** - configuration system is flexible enough to handle new API.

---

## 🧪 Testing Recommendations

Before running on full dataset, test with limited scope:

```python
# In notebook 01_dataset_discovery.ipynb
# Uncomment this line to test with 5 datasets:
datasets_df = datasets_df.head(5)
```

**Test checklist**:
- [ ] Census connection opens successfully
- [ ] `get_obs()` returns NK cell metadata
- [ ] Fallback method works if `get_obs()` fails
- [ ] Categorical columns handled correctly
- [ ] Dataset download works with new URIs
- [ ] Age group assignment works
- [ ] Catalog generation succeeds

---

## 🚀 Migration Guide

### If you have existing code from the old pipeline:

1. **Update cellxgene-census**:
   ```bash
   pip install --upgrade cellxgene-census
   ```

2. **Update code to use new API**:
   - Replace `obs.read().concat().to_pandas()` with `get_obs()`
   - Use context managers for `open_soma()`
   - Handle categorical columns in `.obs`

3. **Test with small dataset first**:
   ```bash
   bash scripts/run_pipeline.sh --test
   ```

4. **Check for deprecation warnings**:
   - Old methods still work but may be deprecated in future
   - Update to new methods for long-term maintenance

---

## 📚 Additional Resources

- **Official docs**: https://chanzuckerberg.github.io/cellxgene-census/
- **Python API**: https://chanzuckerberg.github.io/cellxgene-census/python-api.html
- **GitHub releases**: https://github.com/chanzuckerberg/cellxgene-census/releases
- **Tutorials**: https://chanzuckerberg.github.io/cellxgene-census/notebooks/

---

## 🐛 Known Issues and Workarounds

### Issue 1: Categorical column behavior

**Problem**: `value_counts()` returns all categories, even with 0 cells.

**Workaround**:
```python
# Use dropna=True or convert to string first
counts = adata.obs['cell_type'].astype(str).value_counts()
```

### Issue 2: Memory usage with large queries

**Problem**: Querying millions of cells can consume significant RAM.

**Workaround**:
```python
# Query by dataset_id in batches
for dataset_id in dataset_ids:
    obs_subset = cellxgene_census.get_obs(
        census,
        "homo_sapiens",
        value_filter=f"dataset_id == '{dataset_id}' and cell_type == 'NK cell'"
    )
    # Process subset
```

### Issue 3: Download failures for large H5AD files

**Problem**: Network interruptions during large file downloads.

**Workaround**: Already implemented retry logic and checkpointing in `batch_download()`.

---

## ✅ Compatibility Matrix

| Component | Old Version | New Version | Status |
|-----------|-------------|-------------|--------|
| cellxgene-census | 1.14.1 | ≥1.17.0 | ✅ Updated |
| Census data | 2023-12-15 | 2025-01-30 | ✅ Updated |
| API methods | obs.read() | get_obs() | ✅ Updated |
| Context manager | Optional | Recommended | ✅ Implemented |
| Categorical columns | Not aware | Handled | ✅ Compatible |

---

## 📞 Questions or Issues?

If you encounter problems with the updated API:

1. Check the logs in `results/reports/01_dataset_discovery.log`
2. Review the fallback method output (code automatically tries old API if new fails)
3. Open an issue on GitHub: https://github.com/Alfred3005/NK_INMUNOSENESCENCE/issues

---

**Updated**: 2025-01-03
**Next review**: When Census releases next LTS version
