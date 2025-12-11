---
description: Search the Skills Registry for relevant learnings before starting new work
allowed-tools: Glob, Grep, Read
---

# Skills Registry Advisor

Search the Skills Registry for relevant learnings before starting new work.

## Instructions

1. **Understand the user's goal** from the conversation context: $ARGUMENTS

2. **Search the Skills Registry** at `Skills_Registry/plugins/` for relevant skills:
   - Read all `plugin.json` files to find skills with matching trigger conditions in their descriptions
   - Read corresponding `SKILL.md` files for detailed learnings
   - Look for keywords related to: $ARGUMENTS

3. **Report findings** in this format:

   ### Relevant Skills Found

   For each matching skill:
   - **Skill Name**: [name from plugin.json]
   - **Trigger Match**: [why this skill is relevant]
   - **What Worked**: [key successes from SKILL.md]
   - **What Failed**: [failed attempts table - this is critical]
   - **Recommended Parameters**: [final parameters section]

   ### Recommendations
   Based on these skills, here's what you should consider:
   - [actionable recommendations]

4. **If no relevant skills found**:
   - Inform the user that no matching skills were found
   - Suggest creating a new skill after completing their task using `/retrospective`

## Search Paths
- `Skills_Registry/plugins/training/*/`
- Look for `.claude-plugin/plugin.json` and `skills/*/SKILL.md`
