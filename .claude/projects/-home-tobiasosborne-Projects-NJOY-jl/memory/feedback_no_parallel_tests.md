---
name: Never run parallel Julia test processes
description: Multiple Julia processes competing for precompilation cache causes silent failures
type: feedback
---

NEVER run multiple Julia test/precompile processes simultaneously in the same project.

**Why:** Julia's precompilation cache (~/.julia/compiled/) is not safe for concurrent access. Multiple `Pkg.test()` or `julia test/runtests.jl` processes corrupt each other's cache, causing cryptic "failed to precompile" errors that look like real test failures but aren't. This wasted significant time in the NJOY.jl session.

**How to apply:**
- When launching parallel agents, ensure at most ONE agent runs Julia tests
- Other agents should only read/write files, not execute Julia
- After all agents complete, run exactly ONE clean test: `rm -rf ~/.julia/compiled/v1.12/NJOY* && julia --project=. -e 'using Pkg; Pkg.test()'`
- Never background multiple test runs
