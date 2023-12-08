import argparse
import json
import logging
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from kbase.auth import check_token
from kbase.workspace import Workspace, WorkspaceObjectId
from pathlib import Path
from typing import Any, Tuple

KB_AUTH_TOKEN = "KB_AUTH_TOKEN"
LOG_FILE_PATH = Path(__file__).parent / "genome_downloader.log"

class DownloaderConfig:
    auth_endpoint: str
    service_base_endpoint: str
    workspace_endpoint: str
    genome_types: list[str]

    def __init__(self, config: dict[str, Any]):
        self.auth_endpoint = config["auth_token_endpoint"]
        self.service_base_endpoint = config["service_endpoint"]
        self.workspace_endpoint = self.service_base_endpoint + "/ws"
        self.genome_types = config["genome_type"]

    @classmethod
    def from_file(cls, config_file: Path) -> "DownloaderConfig":
        with open(config_file) as config_in:
            return cls(json.load(config_in))

def load_args() -> dict[str, Any]:
    """
    Loads up the command line arguments.
    Returns a dictionary with:
    ws_id: the workspace id to download genomes from
    out_dir: the output directory for the downloads
    restart: if True, should restart from previous downloads
    """

    desc = """
        This program downloads all genome objects from a KBase narrative.
        Currently, it only downloads the protein feature subset of the genomes.
        Note that the KB_AUTH_TOKEN environment variable should be set with a valid
        auth token for downloading from a private narrative.

        It makes a set of directories, one for each genome, under the given output
        directory. If the directory doesn't exist, this will raise an error. Any
        existing genome files there might be overwritten.

        e.g.
        simple_genome_downloader -n 12345 -o /home/me/my_12345_genomes
    """
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-n', dest='narrative', type=str, required=True, help="Narrative id (for example, 49058 from https://narrative.kbase.us/narrative/49058)")
    parser.add_argument('-o', dest='out_directory', type=str, required=True, help="Output directory")
    parser.add_argument('--restart', dest='restart', action="store_true", required=False, help="restart a partly complete or failed run")
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    out_dir = Path(args.out_directory)
    if not out_dir.is_dir():
        raise ValueError(f"{out_dir} must be a valid directory")
    return {
        "ws_id": args.narrative,
        "out_dir": out_dir,
        "restart": args.restart
    }

def load_config() -> DownloaderConfig:
    """
    Loads and returns the config file. Expects config.json to be present in the same directory
    as the downloader script.
    """
    config_file_path = Path(__file__).parent / "config.json"
    if not config_file_path.is_file():
        raise RuntimeError(f"Missing config file {config_file_path}")
    return DownloaderConfig.from_file(config_file_path)

def get_genome_list(ws: Workspace, ws_id: int, genome_types: list[str], out_file_path: Path) -> list[WorkspaceObjectId]:
    """
    Gets the UPAs for all genomes in the given workspace and dumps them to the given
    file as text, one UPA per line.

    :param ws: the Workspace object for downloading
    :param ws_id: the workspace id to scan for genomes
    :param genome_types: the KBase genome typed object names to use for downloading
    :param out_file_path: the full length file path to the text file for writing
    """
    # get the list
    # dump to file, one per line
    # return the list
    genomes_list = []
    for genome_type in genome_types:
        genomes_list += ws.get_object_upas(ws_id, genome_type)
    with open(out_file_path, "w") as outfile:
        outfile.writelines([f"{upa}\n" for upa in genomes_list])
    return genomes_list

def write_genome_info(genome_data: dict, genome_path: Path) -> None:
    """
    Writes out the genome object info structure to a JSON file.
    :param genome_data: the full genome object from the Workspace, expected to have an "info" key.
    :param genome_path: the path to the output JSON file.
    """
    genome_info = Workspace.obj_info_to_json(genome_data["info"])
    with open(genome_path / "info.json", "w") as info_json:
        json.dump(genome_info, info_json, indent=4)

def build_feature_name(feat) -> str:
    """
    Builds the name of the given feature as a single line for a FASTA header.
    Uses the following structure of whichever keys are available:
    functions=fun1,fun2; functional_descriptions=fdesc1,fdesc2; aliases=ali1,ali2; db_xrefs=db1:ref1,db2:ref2

    This looks for the keys "functions", "functional_descriptions", "aliases", and "db_xrefs", any of which
    may or may not be present.
    """
    name_parts = []
    if feat.get("functions"):
        name_parts.append(f"functions={','.join(feat['functions'])}")
    if feat.get("functional_descriptions"):
        name_parts.append(f"functional_descriptions={','.join(feat['functional_descriptions'])}")
    if feat.get("aliases"):
        name_parts.append(f"aliases={','.join([a[1] for a in feat['aliases']])}")
    if feat.get("db_xrefs"):
        db_xref_info = (f"{db_xref[0]}:{db_xref[1]}" for db_xref in feat["db_xrefs"])
        name_parts.append(f"db_xrefs={','.join(db_xref_info)}")
    return "; ".join(name_parts)

def construct_protein_list(feature_list: list[dict]) -> Tuple[list[SeqRecord], list[str]]:
    """
    Constructs a list of BioPython SeqRecords, and a list of feature ids with no available
    protein sequence.
    The SeqRecords are given the id of the feature, and a name constructed based on
    feature information (see build_feature_name above).
    """
    protein_list = []
    missing = []
    for feat in feature_list:
        if not feat.get("protein_translation"):
            missing.append(feat.get("id", "unknown_feature_id"))
            continue
        protein_list.append(SeqRecord(
            Seq(feat["protein_translation"]),
            id=feat.get("id", "unknown_feature_id"),
            description=build_feature_name(feat)
        ))
    return (protein_list, missing)

def write_genome_protein_fasta(genome_data: dict, genome_path: Path) -> None:
    """
    Writes the given genome object data out to two files in genome_path.
    Up to two files will be written.
    1. the Genome file, with the name of the genome object .faa, as a protein fasta file.
    2. a text file "missing_features_<genome_name>.txt" which is a list of feature
        ids in the genome without any protein sequence.
    :param genome_data: the Genome data object from the Workspace. Expected to have
        a "data" key that holds all the data, and an "info" key with metadata
    :param genome_path: the Path to the directory where files will be written.
    """
    if "features" in genome_data["data"]:
        (protein_seqs, missing) = construct_protein_list(genome_data["data"]["features"])
    elif "cdss" in genome_data["data"]:
        (protein_seqs, missing) = construct_protein_list(genome_data["data"]["cdss"])
    else:
        raise ValueError("feature information not found")
    file_name = f"{genome_data['info'][1]}.faa"
    with open(genome_path / file_name, "w") as genome_out:
        SeqIO.write(protein_seqs, genome_out, "fasta")
    if len(missing):
        with open(genome_path / ("missing_features_" + file_name), "w") as missing_out:
            missing_out.write("\n".join(missing))
        logging.warning(f"genome had no protein translation for {len(missing)} features")

def download_genomes(ws: Workspace, genome_list: list[WorkspaceObjectId], out_file_path: Path, format: str="faa"):
    """
    Downloads genomes from the KBase workspace and dumps them to file.
    The only format working right now is faa - Protein FASTA file.
    """
    if not format == "faa":
        ValueError("Only downloading in FASTA amino acid format right now")
    data_paths = [
        "dna_size",
        "cdss/[*]/id",
        "cdss/[*]/functions",
        "cdss/[*]/functional_descriptions",
        "cdss/[*]/aliases",
        "cdss/[*]/db_xrefs",
        "cdss/[*]/protein_translation",
        "features/[*]/id",
        "features/[*]/functions",
        "features/[*]/functional_descriptions",
        "features/[*]/aliases",
        "features/[*]/db_xrefs",
        "features/[*]/protein_translation",
    ]
    num_genomes = len(genome_list)
    for idx, upa in enumerate(genome_list):
        cur_genome = idx + 1
        logging.info(f"{cur_genome}/{num_genomes} fetching genome {upa} ")
        try:
            genome_data = ws.get_objects([upa.upa], data_paths)[0]
            logging.info(f"{cur_genome}/{num_genomes} got genome {upa}")
            genome_path = out_file_path / str(upa.ws_id) / str(upa.obj_id)
            os.makedirs(genome_path, exist_ok=True)
            logging.info(f"{cur_genome}/{num_genomes} writing {upa} to disk")
            write_genome_info(genome_data, genome_path)
            write_genome_protein_fasta(genome_data, genome_path)
            logging.info(f"{cur_genome}/{num_genomes} done saving {upa}")
        except Exception as err:
            err_str = f"{cur_genome}/{num_genomes} Unable to fetch genome {upa}! {err}"
            logging.error(err_str)
            print(err_str, file=sys.stderr)

def get_restarted_genome_list(genomes_file_path: Path, out_dir: Path, format) -> list[WorkspaceObjectId]:
    """
    This works as a checkpoint manager for the genome downloader.
    On the first pass, a genomes list file is created - the list of all genome UPAs to download.
    If there's a failure of some sort, this is expected to get invoked (ssee the --restart flag).
    It scans that list of genome UPAs, then looks through the directories to see which ones were
    downloaded successfully.
    It then constructs a list of remaining genomes to download, plus the last one that was
    already downloaded, in case it was only a partial download.

    :param genomes_file_path: the Path to the genomes list file, expects to see a single UPA per line.
    :param out_dir: the directory where downloaded genomes are stored. In this directory, there should be
        another dir for the workspace id, and one for each genome object id.
    :param format: the format for genomes being downloaded (should be faa for now), used to
        look for file extensions.
    """
    initial_genomes_list = []
    genomes_list = []
    with open(genomes_file_path) as genomes_file:
        for line in genomes_file:
            initial_genomes_list.append(line.strip())
    restarting_msg = f"Restarting. Found {len(initial_genomes_list)} genomes total."
    print(restarting_msg)
    logging.info(restarting_msg)
    last_done_idx = 0
    for idx, genome in enumerate(initial_genomes_list):
        upa = WorkspaceObjectId.from_upa(genome)
        genome_dir = out_dir / str(upa.ws_id) / str(upa.obj_id)
        is_done = False
        if genome_dir.is_dir() and (genome_dir / "info.json").is_file():
            # check for genome file
            for handle in os.scandir(genome_dir):
                if handle.name.endswith(f".{format}"):
                    is_done = True
                    last_done_idx = idx
        if not is_done:
            genomes_list.append(genome)
    # in case the last one downloaded failed, use it, too.
    genomes_list.append(initial_genomes_list[last_done_idx])
    print(f"{len(initial_genomes_list) - len(genomes_list)} genomes appear to be complete.")
    return [WorkspaceObjectId.from_upa(upa) for upa in genomes_list]

def run_genome_downloader(config: DownloaderConfig, args: dict[str, Any], format="faa") -> None:
    """
    Run the genome downloader!
    1. check that the auth token in KB_AUTH_TOKEN is valid.
    2. build a genomes list, either from scratch or from restarting.
    3. Pull the genomes in that list, one at a time (for now, unless that proves to take too long).
    """
    if not os.environ[KB_AUTH_TOKEN]:
        raise RuntimeError("The KB_AUTH_TOKEN environment variable must be set")

    token = os.environ[KB_AUTH_TOKEN]
    check_token(token, config.auth_endpoint)

    ws_id = args["ws_id"]
    ws = Workspace(token, config.workspace_endpoint)
    genomes_list_file = args["out_dir"] / "genomes_list.txt"
    genomes_list = None
    if args["restart"] and genomes_list_file.is_file():
        # restart from here. get amended genomes list
        genomes_list = get_restarted_genome_list(genomes_list_file, args["out_dir"], format)
    elif args["restart"] and not genomes_list_file.is_file():
        not_restarting_msg = "Restart flag present, but no genomes list file found, starting from scratch..."
        print(not_restarting_msg)
        logging.info(not_restarting_msg)
    if genomes_list is None:
        genomes_list = get_genome_list(ws, ws_id, config.genome_types, genomes_list_file)
    print(f"downloading {len(genomes_list)} genomes")
    print(f"See logfile {LOG_FILE_PATH} for details")
    download_genomes(ws, genomes_list, args["out_dir"], format=format)

if __name__ == "__main__":
    logging.basicConfig(
        filename=LOG_FILE_PATH,
        format="%(asctime)s %(message)s",
        encoding="utf-8",
        level=logging.DEBUG
    )
    config = load_config()
    args = load_args()
    run_genome_downloader(config, args)
