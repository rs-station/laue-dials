# Get package versioning information


def split_stills_by_image(expts, refls):
    """
    Splits the data into lists ordered by image.

    This function takes the ExperimentList and reflection_table provided
    and returns two lists. The first list contains ExperimentList objects
    for each image in the order of existence in the input object, and the
    second list contains corresponding reflection_table objects, including
    any unindexed reflections.

    Args:
        expts (dxtbx.model.ExperimentList) : Input ExperimentList to split by image.
        refls (dials.array_family.flex.reflection_table : Input reflection_table to split by image.

    Returns:
        new_expts (dxtbx.model.ExperimentList) : Output ExperimentList split by image.
        new_refls (dials.array_family.flex.reflection_table : Output reflection_table split by image.
    """
    new_expts = []
    new_refls = []

    unindexed_refls = refls.select(refls["id"] == -1)
    for i in range(len(expts)):
        image_expts = expts[i]
        image_refls = refls.select(refls["id"] == i)
        image_refls.extend(unindexed_refls.select(unindexed_refls["imageset_id"] == image_refls["imageset_id"][0]))

        new_expts.append(image_expts)
        new_refls.append(image_refls)

    return new_expts, new_refls


def laue_version():
    """
    Print the versions of DIALS and laue-dials packages.

    This function retrieves the versions of the DIALS and laue-dials packages and prints
    them to the terminal.

    Args:
        None

    Returns:
        None
    """
    try:
        width = os.get_terminal_size().columns
    except OSError:
        width = 65
    dials_version = version("dials")
    laue_dials_version = version("laue-dials")

    print("-" * width)
    print(f"DIALS version " + dials_version + ".")
    print(f"laue-dials version " + laue_dials_version + ".")
    print("-" * width)
