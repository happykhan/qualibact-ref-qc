import requests

def download_from_google_drive(file_id, destination):
    URL = "https://drive.google.com/uc?export=download"
    session = requests.Session()

    response = session.get(URL, params={'id': file_id}, stream=True)
    token = get_confirm_token(response)

    if token:
        response = session.get(
            URL,
            params={'id': file_id, 'confirm': token},
            stream=True
        )

    save_response_content(response, destination)


def get_confirm_token(response):
    for key, value in response.cookies.items():
        if key.startswith("download_warning"):
            return value
    return None


def save_response_content(response, destination):
    chunk_size = 32768
    with open(destination, "wb") as f:
        for chunk in response.iter_content(chunk_size):
            if chunk:
                f.write(chunk)


file_id = "1lq5zFmru0_xPXgBxBJ5IxZFboiLha6AI"
destination = "got.tar"
download_from_google_drive(file_id, destination)
