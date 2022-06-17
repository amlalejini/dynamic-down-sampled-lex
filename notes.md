# Notes

Ryan's min-max implementation: <https://github.com/ryanboldi/propeller/blob/49d12873350c39cd288e6300bd6fa3af5512f1e7/src/propeller/downsample.cljc#L45-L65>

## Todos
### Selection analyzer

- [ ] Web tool?
- [ ] Add csv parser to submodules

### Diagnostics
- [ ] add output resolution
- [ ] persistent score matrix for population => stop copying repeatedly for every selection event
- [ ] track distribution of test case ids that come up first?
- [ ] pin empirical branch
- [ ] implement toy maxmin sampling in python to play with
- [ ] implement data aggregation script
- [ ] add run config output file
- [x] Stick diagnostics into their own namespace
- [x] add vanilla lexicase, use if epsilon = 0
- [x] add down-sampled lexicase
- [x] implement min-max dynamic down-sampling
- [x] add down-sampled epsilon lexicase
- [x] test maxmin sampling
- [x] implement fair maxmin sampling (i.e., eval % of parents on all test cases, use that info)
- [x] implement lexicase where lead tests are evenly distributed
- [x] switch diagnostic selection to string parameter
- [x] build out HPC submission scripts