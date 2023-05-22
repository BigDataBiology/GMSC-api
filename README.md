# GMSC API

## API Endpoints

- `https://{{base_url}}/v1/seq-info/{{gmsc_id}}`

Where {{gmsc_id}} is of the form `GMSC10.100AA.xxx_xxx_xxx` or `GMSC10.90AA.xxx_xxx_xxx`.

Returns

```json
{
    "id": "GMSC10.xxAA.xxx_xxx_xxxx",
    "nucleotide": "ATC...",
    "aminoacid": "MAV...",
    "taxonomy": "s__Bacteroides_vulgatus",
    "habitat": "human gut"
}
```

### Sequence search interface (non-public)

- `https://{{base_url}}/v1/seq-search`

POST:

- `seq`: `"MVAAK..."`

Returns

```json
{
    "status": "message (normally 'Ok')",
    "search-id": "xxxxx"
}
```

- `https://{{base_url}}/v1/seq-search-results/{{search_id}}`

Returns

```json
{
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

