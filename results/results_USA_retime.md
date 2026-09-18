# Clean re-timing (one query at a time, warm cache)

## ev

| Method | Mean query (s) | Mean nodes expanded | Correct |
|---|---:|---:|---:|
| johnson | 11.890 | 12,832,909 | 10/10 |
| domain | 10.938 | 9,545,337 | 10/10 |
| alt | 3.666 | 1,330,732 | 10/10 |
| max | 4.138 | 1,329,583 | 10/10 |
| bidir-max | 2.812 | 499,533 | 10/10 |
| naive | 12.121 | 10,421,359 | 5/10 |
| 2 cores | 1.171 | 514,677 | 10/10 |

## shifted

| Method | Mean query (s) | Mean nodes expanded | Correct |
|---|---:|---:|---:|
| johnson | 9.250 | 12,670,325 | 10/10 |
| alt | 1.876 | 916,938 | 10/10 |
| bidir-alt | 1.546 | 410,556 | 10/10 |
| 2 cores | 0.829 | 428,622 | 10/10 |
