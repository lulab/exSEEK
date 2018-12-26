# Cluster configuration

**cluster configuration**: `config/cluster.yaml`

Here is an example configuration:

```yaml
__default__:
  queue: queue
  name: {rule}.{wildcards}
  stderr: logs/cluster/{rule}/{wildcards}.stderr
  stdout: logs/cluster/{rule}/{wildcards}.stdout
  threads: {threads}
  resources: span[hosts=1]
```

**cluster command**: `config/cluster_command.txt`

```
bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}
```

**Commonly used parameters**

| Parameter | Description |
| ------ | ----------- |
| `__default__` | Rule name (`__default__`) for default configuration) | 
| `queue` | Queue name (required) |
| `name` | Job name |
| `stderr` | Log file for standard error |
| `stdout` | Log file for standard output |
| `threads` | Number of parallel threads for a job |
| `resources` | Resource requirements. `span[hosts=1]` prevents parallel jobs from being submitted to different nodes |

Refer to the snakemake [documentation](https://snakemake.readthedocs.io/en/stable/).
