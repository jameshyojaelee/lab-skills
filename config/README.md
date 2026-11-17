# Private configuration

Skills that file results somewhere (notebook databases, chat channels, storage paths) read their
destination from a private YAML outside this repo. Point to it with an environment variable:

    export LAB_SKILLS_CONFIG=~/path/to/private/lab-skills-config.yaml

Keep that file out of git (or in a separate private repo). Shape: see `example.yaml`.
Agents: read the file named in `$LAB_SKILLS_CONFIG`; if it is unset, ask the user for the destination.
Never copy values from it into this repo.
