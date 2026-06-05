# Issue 07: Configurable WSClean Command

Type: AFK

## What to Build

Add a WSClean command configuration value and CLI flag. The command should default to `wsclean`, but users should be able to pass a container invocation such as `singularity exec /mnt/software/containers/wsclean-3.10-dysco.sif wsclean`.

## Acceptance Criteria

- [ ] A WSClean command configuration value defaults to `wsclean`.
- [ ] `--wsclean-command` sets the WSClean command from the CLI.
- [ ] The configured command is parsed into argv and executed without shell execution.
- [ ] Execution uses an explicit working directory rather than changing process-global current directory.
- [ ] Tests verify argv construction for both plain `wsclean` and the Singularity invocation.

## Blocked By

- Issue 06: Explicit Imager CLI

## Suggested Commit Message

`feat: configure WSClean command execution`
