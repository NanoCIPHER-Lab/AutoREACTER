![AutoREACTER Logo](_static/logo.png)

## Cleanup Mode

This feature allows you to clean up cached data based on retention time or remove all cached runs.  
It is useful for freeing disk space and maintaining a clean working environment.

Cleanup can be performed via both the CLI and Jupyter Notebook.

---

### Cache cleanup via Jupyter Notebook

In this cell, you just have to specify the number of dates backward or `all` as the mode.

```python
# 'all' → remove all cached data
# int   → keep last N days (e.g., mode=3 keeps last 3 days)
mode = 'all'
```

### Cache cleanup via (CLI)

Run the interactive cleanup utility from the command line:
```python
# Delete runs older than N days (e.g., 7, 30)
python AutoREACTER.py --cleanup N
python AutoREACTER.py -c N
```
```python
# Delete all cached runs
python AutoREACTER.py --cleanup all
python AutoREACTER.py -c all
```

**Notes:**

1. In each mode, the yyyy-mm-dd folder for today will not be deleted
2. Use all with caution — this will permanently delete all cached runs
3. Using a numeric value (e.g., 7) is safer for routine cleanup
4. Jupyter Notebook examples are available in the examples directory
4. If not, you can simply delete manually the `cache/` folder in the root directory.