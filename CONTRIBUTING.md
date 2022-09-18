# Guidelines for contributing

## Issue tracker

Please use the project's GitHub [issue tracker][res-issue-tracker] to:

- find issues to work on
- report bugs
- propose features
- discuss future directions

## Submitting issues

Please choose a template when submitting an issue: choose the [Bug report
template][res-bug-report] only when reporting bugs; for all other issues,
choose the [Feature request template][res-feature-request]. Please follow the
instructions in the template.

You do not need to worry about adding labels or milestones for an issue, the
project maintainers will do that for you. However, it is important that all
issues are written concisely, yet with enough detail and with proper
references (links, screenshots, etc.) to allow other contributors to start
working on them. For bug reports, it is essential that they include
reproducible examples.

## Development

To make it easier for everyone to maintain, read and contribute to the code,
as well as to ensure that the code base is robust and of high quality, we
would kindly ask you to stick to conform to the used code, docstring and commenting style within the project to maintain consistency.
Please stick to the [`black`][res-py-black] Python code formatting style.
For R scripts please use [`styler`][res-r-styler] and for bash scripts: [`beautysh`][res-bash-beautysh].

## Merging your code

Here is a check list that you can follow to make sure that code merges
happen smoothly:

1. [Open an issue](#submitting-issues) first to give other contributors a
   chance to discuss the proposed changes (alternatively: assign yourself
   to one of the existing issues)
2. Fork the repository, clone it & create a feature branch off of the default branch
   (never commit changes to protected branches directly) and implement your
   code changes
3. Ensure that your coding style is in line with the guidelines described above
4. Ensure that all tests configured in the [continuous integration][res-ci-cd] (CI) pipeline pass without
   issues
5. If necessary, clean up excessive commits with `git rebase`; cherry-pick and
   merge commits as you see fit; use concise and descriptive commit messages
6. Push your clean feature branch to the remote; make
   sure the CI pipeline passes
7. Issue a pull request from your fork against the default branch in the original repository; follow the instructions in
   the [template][res-pull-request]; importantly, describe your changes in
   detail, yet with concise language, and do not forget to indicate which
   issue(s) the code changes resolve or refer to; assign a project maintainer
   to review your changes


[res-issue-tracker]: <https://github.com/gruber-sciencelab/MAPP/issues>
[res-bug-report]: .github/ISSUE_TEMPLATE/bug_report.md
[res-feature-request]: .github/ISSUE_TEMPLATE/feature_request.md
[res-py-black]: <https://github.com/psf/black>
[res-r-styler]: <https://github.com/r-lib/styler>
[res-bash-beautysh]: <https://github.com/lovesegfault/beautysh>
[res-ci-cd]: <https://en.wikipedia.org/wiki/Continuous_integration>
[res-pull-request]: PULL_REQUEST_TEMPLATE.md