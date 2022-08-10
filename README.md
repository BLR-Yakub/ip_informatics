# IP_quality_informatics
This project aims at automating the process of sequence matching with the sequences available on publicly available patent databases such as WIPO and USPTO.

#### Environment and Packages required
Python >=3.8
Selenium ```pip install selenium```  
Biopython ```pip install biopython```  
Flask ```pip install Flask```  
Chrome webdriver latest version

#### Usage

* Use the jupyter notebook get_sequences.ipynb to extract the sequences from wipo or uspto database
* Use the jupyter notebook velocity_blast.ipynb to do the velocity blast for a single or multiple sequences
* Web app can be used in the host computer by first unzipping the fetched_uspto_all_final and wipo final files into csv and making three files patent_bridge_app.py, index.html, app_css.css and then running command ```python patent_bridge_app.py``` in the terminal.

#### Issues and limitations
* Package is not built yet, so run the code by installing the necessary environments in the host computer.
* Web App is not hosted yet, so it needs to be run on the local server
