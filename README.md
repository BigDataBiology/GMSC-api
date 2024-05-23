# GMSC API

## API Endpoints

- `https://{{base_url}}/v1/seq-info/{{gmsc_id}}`

Where {{gmsc\_id}} is of the form `GMSC10.100AA.xxx_xxx_xxx` or `GMSC10.90AA.xxx_xxx_xxx`.

Returns

```json
{
    "id": "GMSC10.xxAA.xxx_xxx_xxxx",
    "nucleotide": "ATC...",
    "aminoacid": "MAV...",
    "taxonomy": "s__Bacteroides_vulgatus",
    "habitat": "human gut",
    "quality": {
        "antifam": true,
        "terminal": true,
        "rnacode": 0.9,
        "metat": 1,
        "metap": 1,
        "riboseq": 0.9
    }
}
```

Note that the `quality` field is only present for 90AA sequences.

- `https://{{base_url}}/v1/seq-info-multi/`

This is a `POST`-only endpoint, expecting a JSON package consisting of a
dictonary with an entry `seq_ids`, which is a list of strings (identifiers).
For example:

```json
{
    "seq_ids": [
            "GMSC10.90AA.123_456_789",
            "GMSC10.90AA.123_456_790",
            ...]
}
```

Returns a list of entries like the outputs of `seq-info`.

- `https://{{base_url}}/v1/seq-filter/`

`POST` endpoint, with arguments:

- `hq_only`: boolean. _optional_ (only active for 90AA)
- `habitat`: str. _mandatory_
- `taxonomy`: str. _optional_
- `quality_antifam`: boolean. _optional_
- `quality_terminal`: boolean. _optional_
- `quality_rnacode`: float. _optional_
- `quality_metat`: integer. _optional_
- `quality_metap`: integer. _optional_
- `quality_riboseq`: float. _optional_

`habitat` is treated as a comma separated list (_e.g._, you can use `marine,freshwater` to match all the entities that are present in both marine **and** freshwater).

`taxonomy` is a substring match so you can pass any taxonomic level (_e.g._, passing `o__Pelagibacterales` will match `d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Pelagibacterales;f__Pelagibacteraceae;g__AAA240-E13`).

Returns

```json
{
    "status":"Ok",
    "results": [
            {
                "habitat":"marine,plant associated,sediment",
                "seq_id":"GMSC10.90AA.000_013_322",
                "taxonomy":"d__Bacteria"},
                ....
    ]
}
```

At most 1,001 entries are returned.

- `https://{{base_url}}/v1/cluster-info/{{gmsc_90AA_id}}`

Returns the membership of the given cluster. At most **20 results** are _thick_ (meaning that metadata is also returned). For the rest, only identifiers are returned. Example output

```json
{
	"status":" Ok",
	"cluster": [
		{
			"aminoacid":"MAAAGFLIVSFKPFEKPSRNAATTAGFSAENFEFTMIALPYSLRP",
			"habitat":"soil",
			"nucleotide":"ATGGCCGCGGCCGGATTCTTGATCGTGTCCTTCAAGCCTTTCGAGAAGCCTTCGAGAAACGCCGCGACGACGGCCGGCTTCTCGGCCGAGAATTTCGAGTTCACGATGATCGCGCTGCCGTACAGCTTGAGACCGTAA",
			"seq_id":"GMSC10.100AA.547_444_661",
			"taxonomy":"d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__VAZQ01;s__VAZQ01 sp005883115"
		}, ...
		]
}
```


### Sequence search interface (non-public interface)

**NOTE**. These are not recommended for public use. For large-scale analyses,
we recommend you use the
[GMSC-mapper](https://github.com/BigDataBiology/gmsC-mapper/) command line tool
locally. Public API endpoints will be maintained for the long-term. No such commitment
is made for endpoints marked `internal`. _You have been warned_.

- `https://{{base_url}}/internal/seq-search` (`POST`)

Arguments:

- `sequence_faa`: FASTA formatted set of sequences
- `is_contigs`: bool (when `True`, inputs are assumed to be DNA contigs)

Returns

```json
{
    "status": "message (normally 'Ok')",
    "search-id": "xxxxx"
}
```

- `https://{{base_url}}/internal/seq-search/{{search_id}}`

Returns

```json
{
    "search_id": "str",
    "status": "str",
    "results": [
        {
            "query_id": "query_1",
            "aminoacid": "MHEDVIQFARNEVWSLV....",
            "taxonomy": "s__Bacteroides_vulgatus",
            "habitat": "human gut",
            "hits": [
                { "id": "GMSC10.xxAA.xxx_xxx_xxxx",
                  "e_value": "2.1e-23",
                  "aminoacid": "MHEELIQFARNEV...",
                  "identity": "98.4"
                }, ...
            ]
        }, ...]
```

`status` will be one of `Running` (if the results are not yet ready), `Done`,
or `Expired`. In the case of `Done`, the `results` field will be filled in.

## Install & Testing

Dependencies

- `flask`
- `numpy`
- `pandas`

Running this (in test mode) can be done with

```bash
python -m flask run
```

Testing can be done with `curl`:

```bash
curl http://127.0.0.1:5000/v1/seq-info/GMSC10.100AA.000_000_002
```

These examples assume you are running the test version on
`http://127.0.0.1:5000/`. Adapt as necessary.

Searching requires using `POST` and a FASTA file. For example, if you have the
file `example.faa`, you can use

```bash
curl -X POST --form "sequence_faa=$(cat example.faa)" http://127.0.0.1:5000/internal/seq-search/
```

The output will look something like this:

```json
{"search_id":"1-jmgi","status":"Ok"}
```

You can later use the given ID (in this case `1-jmgi`, but it will be different
every time the app runs) to retrieve the results:

```bash
curl http://127.0.0.1:5000/internal/seq-search/1-jgmi
```

Results will look like one of the following

1. `{"search_id":"1-jmgi","status":"Running"}`
2. `{"search_id":"1-jmgi","status":"Done", results":[...]}`
3. `{"search_id":"1-jmgi","status":"Expired"}`

Search ID are of the form `#-xxxx` where `#` is just an index counting up and
`xxxx` is a random string.


## Indexing

Indexing is done by the `make-indices.py` [Jug](https://jug.readthedocs.io/)
script. It expects FASTA and other files to be present in the `gsmc-db`
subdirectory.

