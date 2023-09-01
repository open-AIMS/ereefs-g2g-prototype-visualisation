# Copyright Eric Lawrey, Australian Institute of Marine Science
# This script download the shapefiles needed to make the basemap in the plots.
# The data can be manually download in a browser from:
# https://nextcloud.eatlas.org.au/s/RGwTFcLtmPApEcQ/download
import urllib.request
import zipfile
import os

def urlretrieve(url, path, reporthook=None, user_agent=None):
    request = urllib.request.Request(url)
    if user_agent:
        request.add_header("User-agent", user_agent)

    with urllib.request.urlopen(request) as response, open(path, 'wb') as out_file:
        meta = response.info()
        file_size = meta.get("Content-Length")
        block_size = 8192*32
        read_so_far = 0
        
        # Convert file_size to int if it's provided, else None
        file_size = int(file_size) if file_size else None

        while True:
            buffer = response.read(block_size)
            if not buffer:
                break

            read_so_far += len(buffer)
            out_file.write(buffer)

            if reporthook is not None:
                # Call reporthook with available data, if no file_size, use read_so_far as a placeholder
                reporthook(read_so_far // block_size, block_size, file_size)

def show_progress(count, block_size, total_size):
    if total_size:
        percent_complete = (count * block_size) * 100 / total_size
        print(f"Downloaded {count * block_size} of {total_size} bytes ({percent_complete:.2f}% complete)")
    else:
        print(f"Downloaded {count * block_size} bytes (total size unknown)")

    
def download_and_unzip(url, destination_folder):
    """
    Downloads and unzips a file from a specified URL using only Python's built-in libraries.
    This version provides feedback during the download and specifies a user agent.

    Parameters:
    - url: The URL of the file to be downloaded.
    - destination_folder: Folder where the downloaded file will be saved and unzipped.
    """
    
    # Ensure the destination folder exists
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)

    # Define a temporary path for the downloaded file
    zip_file_path = os.path.join(destination_folder, "temp_file.zip")

    # Download the file with progress feedback and user agent
    user_agent = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3"
    headers = {"User-Agent": user_agent}
    request = urllib.request.Request(url, headers=headers)

    # Downloading the file and providing feedback
    #urllib.request.urlretrieve(url, zip_file_path, reporthook=show_progress)
    
    urlretrieve(url, zip_file_path, show_progress, user_agent=user_agent)

    # Unzip the file
    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        zip_ref.extractall(destination_folder)

    # Remove the downloaded zip file
    os.remove(zip_file_path)
    print(f"\nFile from {url} has been downloaded and unzipped to {destination_folder}.")

if __name__ == "__main__":
    # Example usage
    URL = "https://nextcloud.eatlas.org.au/s/RGwTFcLtmPApEcQ/download"
    DESTINATION = "src-data"
    
    download_and_unzip(URL, DESTINATION)
