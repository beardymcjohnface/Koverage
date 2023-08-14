rule build_environment:
    input:
        os.path.join(dir["env"], "{env}.{yaml}")
    output:
        temp(touch(dir["temp"], "{env}.{yaml}.done"))
    wildcard_constraints:
        yaml = "yaml|yml"
    conda:
        os.path.join(dir["env"], "{env}.{yaml}")
