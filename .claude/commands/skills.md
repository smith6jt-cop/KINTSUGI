# List Skills Registry

List all skills in the Skills Registry with their descriptions and trigger conditions.

## Instructions

1. **Scan the Skills Registry** at `Skills_Registry/plugins/`

2. **For each skill found**, read the `plugin.json` and extract:
   - Name
   - Description (including trigger conditions)
   - Author
   - Version

3. **Display as a formatted table**:

   | Skill | Description | Author |
   |-------|-------------|--------|
   | skill-name | Trigger conditions... | Author Name |

4. **Show summary statistics**:
   - Total number of skills
   - Skills by category (based on folder structure in plugins/)

## Search Path
`Skills_Registry/plugins/*/*/`
