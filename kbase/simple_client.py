import requests
import uuid
from typing import Any

def make_kbase_jsonrpc_1_call(
        service_endpoint: str,
        service: str,
        method: str,
        params: list[Any],
        token: str=None) -> Any:
    """
    A very simple JSON-RPC 1 request maker for KBase services.

    If a failure happens, it prints the message from an expected error packet, and
    raises the requests.HTTPError.
    """
    call_id = str(uuid.uuid4())
    json_rpc_package = {
        "params": params,
        "method": f"{service}.{method}",
        "version": "1.1",
        "id": call_id
    }
    headers = {}
    if token is not None:
        headers = {"Authorization": token}
    resp = requests.post(service_endpoint, json=json_rpc_package, headers=headers)
    try:
        resp.raise_for_status()
    except requests.HTTPError:
        print(resp.json().get("error").get("message"))
        raise
    return resp.json()["result"]
