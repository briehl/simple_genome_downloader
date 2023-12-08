import requests

def check_token(token: str, auth_endpoint: str) -> dict[str, str]:
    headers = {
        "Authorization": token
    }
    resp = requests.get(auth_endpoint, headers=headers, allow_redirects=True)
    try:
        resp.raise_for_status()
        return resp.json()
    except requests.HTTPError:
        if resp.status_code == 401:
            raise ValueError("Invalid KBase auth token")
        else:
            raise ValueError(resp.text)
