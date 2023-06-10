# Dynamic Mode Decomposition

## Introduction
The following is a research project that I worked on as a Postdoctoral Scholar at the School of Mathematical and Statistical Sciences at Arizona State University. As a researcher in the field of fluid dynamics, I pursued a data driven point of view, using modern machine learning algorithms on the time series data of fluid flow states.

Fluid dynamics is an important field of study for many industries including: Aerospace, Energy (wind, hydroelectric, fossil fuels / oil), Automotive, Environmental, Pharmaceutical, Chemical and Petrochemical, Food and Beverage, etc. In general, this is implemented as a complex time series of data, where one must rely on mathematical theorems and algorithms to lead explainability and predictability of the phenomenon recognized in computational simulation models. 

## Description
The [Dynamic Mode Decomposition (DMD) algorithm](https://en.wikipedia.org/wiki/Dynamic_mode_decomposition) is a dimensionality reduction algorithm appearing in Schmid 2010 (see [References](#references)). Given simple dynamics (linear flows), it is known that one can use the DMD algorithm to estimate the future fluid velocity data. In this project, we extend this method to apply to weakly nonlinear flows. This project has applications to turbulent flow predictions such as weather forecasts. 






## Files
All code is written in Matlab (version R2022b: 9.13.0.2105380).

- [ ] DMDreport.pdf : A file detailing the project and results.
- [ ] MatlabCode\DMDflow.m : The main Matlab file which generates all images, videos, and results. Change various Boolean values from False to True to explore different sections of the code and produce some nice pictures.
- [ ]  MatlabCode\DMD.m : Helper function implementing DMD. Mostly borrowed from Kutz et. al. 2016 (see [References](#references)).
- [ ] MatlabCode\DMDflow.m :  Helper function reading the fluid flow data from the Cases folder.
- [ ] Cases : This folder contains the raw fluid flow data as binary files (bin01XXX). The data was generated in Fortran using a fourth order Runge–Kutta method. The Fortran file can be found at: [Diablo.f](https://www.damtp.cam.ac.uk/user/jrt51/files.html).

## References
- [ ] Schmid, P. (2010). Dynamic mode decomposition of numerical and experimental data. _Journal of Fluid Mechanics,_  _656_, 5-28. doi:[10.1017/S0022112010001217](https://doi.org/10.1017/S0022112010001217)
- [ ] Kutz, N. , Brunton, S. , Brunton, B. , and Proctor, J. (2016). Dynamic Mode Decomposition. _Society for Industrial and Applied Mathematics,_ Philadelphia, PA. doi:[10.1137/1.9781611974508](https://epubs.siam.org/doi/abs/10.1137/1.9781611974508)






## Suggestions for a good README
Every project is different, so consider which of these sections apply to yours. The sections used in the template are suggestions for most open source projects. Also keep in mind that while a README can be too long and detailed, too long is better than too short. If you think your README is too long, consider utilizing another form of documentation rather than cutting out information.

## Name
Choose a self-explaining name for your project.

## Description
Let people know what your project can do specifically. Provide context and add a link to any reference visitors might be unfamiliar with. A list of Features or a Background subsection can also be added here. If there are alternatives to your project, this is a good place to list differentiating factors.

## Badges
On some READMEs, you may see small images that convey metadata, such as whether or not all the tests are passing for the project. You can use Shields to add some to your README. Many services also have instructions for adding a badge.

## Visuals
Depending on what you are making, it can be a good idea to include screenshots or even a video (you'll frequently see GIFs rather than actual videos). Tools like ttygif can help, but check out Asciinema for a more sophisticated method.

## Installation
Within a particular ecosystem, there may be a common way of installing things, such as using Yarn, NuGet, or Homebrew. However, consider the possibility that whoever is reading your README is a novice and would like more guidance. Listing specific steps helps remove ambiguity and gets people to using your project as quickly as possible. If it only runs in a specific context like a particular programming language version or operating system or has dependencies that have to be installed manually, also add a Requirements subsection.

## Usage
Use examples liberally, and show the expected output if you can. It's helpful to have inline the smallest example of usage that you can demonstrate, while providing links to more sophisticated examples if they are too long to reasonably include in the README.

## Support
Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc.

## Roadmap
If you have ideas for releases in the future, it is a good idea to list them in the README.

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.
