#!/opt/tljh/user/envs/physio/bin/python
import os
import json


def publish_notebook(nb):
    with open(nb) as f:
        text = json.loads(f.read())

    ic = 0
    for cell in text["cells"]:
        if "ExecuteTime" in cell["metadata"]:
            del cell["metadata"]["ExecuteTime"]
        if cell["cell_type"] == "code":
            ic += 1
            cell['execution_count'] = ic
    basename = os.path.splitext(nb)[0]
    with open(f"{basename}_.ipynb", "w") as f:
        json.dump(text, f)

    command = f"""
        jupyter nbconvert --to markdown {basename}_.ipynb;
        rm {basename}_.ipynb
        """

    if os.system(command) == 0:
        print(f"{nb} published successfully.")
    else:
        print("Oops, something went off.")
    return 0


if __name__ == "__main__":
    from sys import argv

    notebooks = argv[1:]
    notebooks = [nb for nb in notebooks if nb.endswith(".ipynb")]
    for nb in notebooks:
        publish_notebook(nb)
