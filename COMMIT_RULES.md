# Commit rules

Format: `[Tag] Short summary in the imperative`

| Tag       | Use for                                                          |
| --------- | ---------------------------------------------------------------- |
| `[Feat]`  | a new module or feature (`include/*.h`, `src/main.cpp`)          |
| `[Test]`  | tests, test runner, mutation testing                             |
| `[Fix]`   | bug fix                                                          |
| `[Tool]`  | `Makefile`, `build_all.bat`, scripts in `tools/`, `viz/`         |
| `[Doc]`   | README and other documents                                       |
| `[Chore]` | `.gitignore`, `.gitattributes`, images, `third_party/`, folders  |

- Summary line: imperative mood, no full stop, about 70 characters or fewer.
- Use the description box for the why and any detail.
- One task per commit. Commit small and often.
- Pull before you start work and before you push.
- Stage only the files of that commit. Do not use "commit all".
- Build and run the tests before committing code: `make` (Linux, macOS) or `build_all.bat` (Windows). All tests must pass.
- Never commit `results/`, `.exe` files, `__pycache__/` or `*_tmp/` folders. They are generated and git-ignored.
- Do not add an AI tool as co-author.

Example:

```
[Feat] Add dp.h with the DP seam finder

Bottom-up, memoised and rolling-row versions, plus seam counting.
Adds tests DP1 to DP9 (23 tests pass in total).
```
