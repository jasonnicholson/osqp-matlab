## To Do
[x] Remove reference to `v1.0.0`. We need a better way to manage version number. Do not hard code and docstrings everywhere.
[ ] This needs so much work. It's inconsistent. 
    - [x] Constants are not used well. - Handled by extracting properties as doubles directly in configuration.
    - Options are used improperly and not consistent.
[x] Update README.md structure of project.

## To Do Later
[ ] Compile on different OS's. I don't trust that this make script is going to work across all platforms. Focus on Linux and Windows. 2020a to 2026a. 
[ ] 

## Questions
[x] Are enumerations slow? Does int32 versus double matter? MATLAB prefers double. - Resolved: Converted into flat struct for the run loops.