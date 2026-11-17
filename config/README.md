# Private configuration

Skills that file results somewhere (notebook databases, chat channels, storage paths) read their
destination from a private YAML outside this repo. Agents look for it in this order:

1. the path in `$LAB_SKILLS_CONFIG`
2. `~/.config/lab-skills/config.yaml` (a file or a symlink to wherever you keep it)

Keep that file out of git (or in a separate private repo). Shape: see `example.yaml`.
Agents: read the first of these that exists; if neither does, ask the user for the destination.
Never copy values from it into this repo.
