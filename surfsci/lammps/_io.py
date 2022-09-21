# +
def get_input_commands(file: str) -> list[str]:
    """
    Reads an input file and returns a list of commands (strings).
    Comments, empty lines, and empty spaces are eliminated.
    """
    with open(file, "r") as of:
        commands = []
        for line in of.readlines():
            # eliminate comments
            if "#" in line:
                end = line.index("#")
                line = line[:end]
            # eliminate empty spaces
            words = line.split()
            cmd = " ".join(words)
            if cmd != "":
                # TODO: support multi-line commands
                if "&" in cmd:
                    raise RuntimeError("Multi-line commands are not supported yet!")
                commands.append(cmd)
    return commands


def write_commands(commands: list[str], file: str) -> None:
    cmds = []
    args = []
    for line in commands:
        split = line.split()
        cmds.append(split[0])
        args.append(" ".join(split[1:]))
    just = max([len(w) for w in cmds])
    with open(file, "w") as of:
        for cmd, arg in zip(cmds, args):
            of.write(f"{cmd.ljust(just+4)}{arg}\n")
