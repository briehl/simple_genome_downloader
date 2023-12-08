# simple_genome_downloader
A really simple downloader for KBase genomes.

Get an auth token using secret KBase powers and set it to the `KB_AUTH_TOKEN` environment variable.
Make a directory for your genomes, like `my_genomes` (but not that, because that's silly).

Then, given a Narratiev id as nar_id (should be a number, look at your narrative URL), run:  
`python simple_genome_downloader -n nar_id -o my_genomes`

This will make the following files:
* `my_genomes/genomes_list.txt` - a list of all UPAs referencing a genome in that Narrative.
* `my_genomes/nar_id/genome_id/info.json` - the JSON object info for that genome.
* `my_genomes/nar_id/genome_id/genome_name.faa` - the genome protein FASTA file.
* `my_genomes/nar_id/genome_id/missing_features_genome_name.txt` - a list of feature ids without protein translation information for that genome.

If the download fails, run it again with the `--restart` flag:  
`python simple_genome_downloader -n nar_id -o my_genomes --restart`  
This will restart the download with the last genome downloaded.
