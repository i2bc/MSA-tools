import requests

prefix = "https://files.rcsb.org/download/"
suffix = ".cif"

url = prefix + "8k81" + suffix


def download_file(url, output_file):
    """download a file but does not use stream so not suitable to download big files"""
    response = requests.get(url)
    print(f"status code : {response.status_code}")

    if response.status_code == 200:
        with open(output_file, "wb") as f:
            f.write(response.content)
    else:
        raise ValueError(f"url returned status_code {response.status_code}")


download_file(url, "test.cif")
