---
description: Save learnings from the current session as a new skill in the registry
allowed-tools: Glob, Grep, Read, Write
---

# Skills Registry Retrospective

Save learnings from the current session as a new skill in the registry.

## Instructions

1. **Gather information** about the session:
   - Skill name (kebab-case): $ARGUMENTS or derive from conversation
   - What was the goal?
   - What environment was used (hardware, software versions)?
   - What worked?
   - What failed and why?
   - What are the final parameters/configuration?

2. **Create the skill folder structure** by copying the template:
   ```
   Skills_Registry/plugins/training/<skill-name>/
   ├── .claude-plugin/
   │   └── plugin.json
   ├── skills/<skill-name>/
   │   └── SKILL.md
   ├── references/
   │   └── .gitkeep
   └── scripts/
       └── .gitkeep
   ```

3. **Fill in plugin.json** with:
   ```json
   {
     "name": "<skill-name>",
     "version": "1.0.0",
     "description": "Specific trigger conditions: (1) <scenario>, (2) <scenario>. Verified on <environment>.",
     "author": {
       "name": "<author>"
     },
     "skills": "./skills",
     "repository": "https://github.com/smith6jt-cop/Skills_Registry"
   }
   ```

4. **Fill in SKILL.md** with the full skill documentation following the template structure:
   - Experiment Overview table
   - Context section
   - Verified Workflow (step-by-step with code blocks)
   - **Failed Attempts table (CRITICAL - most valuable section)**
   - Final Parameters (exact, copy-pasteable)
   - Key Insights
   - References

5. **Create a branch and commit**:
   ```bash
   cd Skills_Registry
   git checkout -b skill/<skill-name>
   git add .
   git commit -m "Add skill: <skill-name>"
   ```

6. **Open a PR** to main using `gh pr create`

## Template Location
`Skills_Registry/templates/experiment-skill-template/`

## Important Rules
- The Failed Attempts table is the most valuable section - be thorough
- Include exact parameters, not vague advice
- Add specific trigger conditions in the description so /advise can find it
- Document the environment for reproducibility
