from karabo.simulation.telescope import Telescope
from karabo.simulation.telescope import OSKARTelescopesWithVersionType, OSKARTelescopesWithoutVersionType
from karabo.simulation import telescope_versions
from karabo.simulation.visibility import *
from thefuzz import process, fuzz



def letters_and_digits(s):
    """Return a string containing only letters and digits from the input string."""
    return ''.join(filter(str.isalnum, s)).lower().replace('versions', '')

def only_letters(s):
    """Return a string containing only letters from the input string."""
    return ''.join(filter(str.isalpha, s)).lower().replace('versions', '')


def get_telescope_version(telescope_name):
    """
    Get the version of the telescope from the name.
    If the name contains a version, return it, otherwise return None.
    """
    # Check if the telescope name contains a version
    version_telescope = None

    if telescope_name in get_args(OSKARTelescopesWithVersionType):
        from karabo.simulation import telescope_versions
        available_versions = []
        available_objects = []
        for obj in dir(telescope_versions):
            # Check if the object is a class

            if isinstance(getattr(telescope_versions, obj), type):
                available_objects.append(getattr(telescope_versions, obj))
                available_versions.append(letters_and_digits(getattr(telescope_versions, obj).__name__))
        myversions = []
        similarity = process.extract(only_letters(telescope_name), available_versions, scorer=fuzz.ratio)
        for version, score in similarity:
            if score > 90:
                # get the index of the version in available_versions
                idx = available_versions.index(version)
                myversions.append(available_objects[idx])

        available_options = []
        available_parents = []
        for version_option in myversions:
            for option in version_option:
                available_parents.append(version_option)
                available_options.append(option)

        if len(available_options) == 0:
            print(f"No versions found for {telescope_name}. Exiting.")
            return None
        else:
            return available_options
        # if len(available_options) == 1:
        #     version_telescope = available_options[0]
        #     print(f"Using version {version_telescope} for telescope {telescope_name}")
        # elif len(myversions) > 1:
            
        #     option = None
        #     for idx, option in enumerate(available_options):
        #         print(f"{idx+1}. {option}")
        #     choice = input("Enter the number of the version you want to use: ")
        #     try:
        #         choice = int(choice)
        #         if choice < 1 or choice > len(myversions):
        #             raise ValueError("Invalid choice")
        #         version_telescope = available_options[choice-1]
        #     except ValueError as e:
        #         print(f"Invalid choice: {e}. Exiting.")
        #         sys.exit(1)

        # return version_telescope




# telescope = Telescope.constructor("EXAMPLE")
# telescope.plot_telescope(file="example_telescope.png")


# available_versions = []
# available_objects = []
# for obj in dir(telescope_versions):
#     # Check if the object is a class

#     if isinstance(getattr(telescope_versions, obj), type):
#         available_objects.append(getattr(telescope_versions, obj))
#         available_versions.append((getattr(telescope_versions, obj).__name__))

# print (available_versions)

#telescope_types = get_args(OSKARTelescopesWithVersionType) + get_args(OSKARTelescopesWithoutVersionType)
telescope_types = get_args(OSKARTelescopesWithoutVersionType)
telescope_versions = get_args(OSKARTelescopesWithVersionType)

for telescope_type in telescope_types:
    if ("SKA" in telescope_type and "MID" in telescope_type):
        print(f"Processing telescope type: {telescope_type}")
        tel = Telescope.constructor(telescope_type)
        tel.plot_telescope(file=f"{telescope_type}.png")
for telescope_name in telescope_versions:
    if "SKA" in telescope_name.upper() and "MID" in telescope_name.upper():
        versions = get_telescope_version(telescope_name)
        if versions is None:
            print(f"No version found for telescope {telescope_name}. Skipping.")
            continue

        for version in versions:
            try:

                # You can instantiate the telescope with the version here if needed

                tel = Telescope.constructor(telescope_name, version=version)
                tel.plot_telescope(file=f"{telescope_type}_v{version}.png")
                print(f"Found version {version} for telescope {telescope_type}")
            except:
                pass


        # You can instantiate the telescope here if needed
        # tel = Telescope.constructor(telescope_type.name, version=telescope_type.version)
        # tel.plot_telescope(file=f"{telescope_type.name}_v{telescope_type.version}.png")


# for version in available_versions:
#     if ("SKA" in version and "Mid" in version):
#         print (f"Found SKA Mid version: {version}")
#         telescope = Telescope.constructor(args.telescope, backend=backend, version=version_telescope)


    # tel = Telescope.constructor(telescope.name, version=telescope.version)
    # tel.plot_telescope(file=f"{telescope.name}_v{telescope.version}.png")

# for telescope in OSKARTelescopesWithVersionType:
#     print(f"Processing {telescope.name} with version {telescope.version.__name__}")
#     # tel = Telescope.constructor(telescope.name, version=telescope.version)
#     # tel.plot_telescope(file=f"{telescope.name}_v{telescope.version}.png")

# for telescope in OSKARTelescopesWithoutVersionType:
#     print(f"Processing {telescope.name} without version")
#     # tel = Telescope.constructor(telescope.name)
#     # tel.plot_telescope(file=f"{telescope.name}.png")