---
name: Repeat rules after every wave
description: Repeat the 3+1 workflow rules after completing each wave to maintain LLM attention on them
type: feedback
---

After completing each wave, explicitly restate the 3+1 agent workflow rules before starting the next wave.

**Why:** LLM attention mechanisms lose focus on instructions that appeared early in context. Repeating the rules keeps them salient and prevents drift into shortcuts.

**How to apply:** After every wave completion, output a block restating: (1) dual proposals required, (2) independent agents, (3) skeptical reviewer with PASS/FAIL, (4) decision log in reports/decisions/, (5) no merge without reviewer PASS.
