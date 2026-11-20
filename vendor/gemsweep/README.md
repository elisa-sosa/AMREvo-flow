# Pinned GEMSWEEP checkout

Run `bin/bootstrap_gemsweep.sh` to create the ignored `checkout/` from the
upstream MIT-licensed GEMSWEEP repository. The helper checks out commit
`d06885f`, initializes upstream submodules, verifies `LICENSE.md`, and applies
the tracked Nextflow-25 compatibility patch. Upstream code and licensing are
not copied into AMREvoFlow or silently updated.
