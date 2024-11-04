
import requests
from bs4 import BeautifulSoup
import os
# Specify the URL of the website you want to scrape
url = 'https://web.engr.oregonstate.edu/~mjb/cs575/'

# Send an HTTP request to the website and get its content
response = requests.get(url)
content = response.text

# Parse the content using BeautifulSoup
soup = BeautifulSoup(content, 'html.parser')

# Find all the links in the website
all_links = soup.find_all('a')

# Filter the links based on some criteria (e.g., containing a specific word)
word_to_filter = 'pp.dpf'
filtered_links = [link['href'] for link in all_links if word_to_filter in link['href']]

# Save the filtered links to a .txt file
output_file = 'filtered_links.txt'
with open(output_file, 'w') as f:
    for link in filtered_links:
        f.write(link + '\n')

print(f'Filtered links saved to {output_file}')

def download_file(url, local_filename):
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return local_filename

# Replace 'your_directory' with the desired directory where you want to save the files
target_directory = 'D:\\Courses\\23Spring\\CS575\\pdfs'

with open('filtered_links.txt', 'r') as file:
    urls = file.readlines()

for url in urls:
    url = url.strip()
    local_filename = url.split('/')[-1]
    
    # Join the target directory with the local filename
    local_filepath = os.path.join(target_directory, local_filename)
    
    print(f'Downloading {url} to {local_filepath}')
    download_file(url, local_filepath)
    print(f'Successfully downloaded {local_filepath}')


