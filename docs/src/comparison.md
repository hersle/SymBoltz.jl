# Comparison with other codes

TODO

| Feature                  | [SymBoltz.jl](https://github.com/hersle/SymBoltz.jl) | [CAMB](https://camb.info/)    | [CLASS](https://lesgourg.github.io/class_public/class.html) |
| :----------------------: | :--------------------------------------------------: | :---------------------------: | :---------------------------------------------------------: |
| **Language**             | Julia                                                | Fortran + Python              | C + Python                                                  |
| **Modularity**           | Every component is fully modular                     | Class-based for simple models | Requires raw code modification in many places               |
| **Approximations**       | None                                                 | Mandatory (?)                 | Mandatory (?)                                               |
| **Speed**                | Fast                                                 | Faster                        | Fastest                                                     |
| **Sensitivity analysis** | Automatic differentiation                            | Finite differences            | Finite differences                                          |