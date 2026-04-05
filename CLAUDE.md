# laue-dials

At the start of every session, read the full contents of
`docs/laue_dials_llm_programming_guidelines.rst`.

## Git Commits

Always pass commit messages as an explicit `-m` string so the full
message is visible before the user approves the tool call:

```
git commit -m "Short subject line

Body text if needed.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

Never use heredoc (`<<'EOF'`) syntax for commit messages.
