---
name: Follow 3+1 agent workflow strictly
description: User expects strict adherence to the dual-proposal + reviewer process from prompt.md — never skip it
type: feedback
---

Always follow the 3+1 agent workflow from prompt.md for every major decision. This is non-negotiable.

**Why:** The user designed a specific multi-agent workflow (Researcher → 2 independent Proposers → Skeptical Reviewer) to ensure quality and catch errors. Skipping steps undermines the whole process.

**How to apply:**
- For every module port, architectural decision, or implementation task listed in the prompt, spawn TWO independent proposer agents
- After comparing and selecting, spawn a skeptical reviewer agent
- Log every decision in reports/decisions/NNN-topic.md
- Never commit code without a reviewer PASS
- If a proposal method fails (e.g., worktree not available), find another way to get an independent second proposal — don't just go with the single one
