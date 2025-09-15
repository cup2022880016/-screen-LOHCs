Screening LOHCs
======
Function introduction
---
  For screening qualified liquid organic hydrogen carriers from PubChem database.
  
Required Python package
---
  Pandas, RDKit, pubchempy, and alfabet.

Use examples
---
  * Open screen_lohc_main.py in the compiler and run it. Enter the required number of carbon atoms and whether heteroatoms and elements are required according to the prompts.
  * Example is the output example after running (the number of carbon atoms input is 8, and there are no heteroatoms). Data contains all the output results.
  * If you need data for each stage, use the to_excel function under each variable in screen_lohc_main.py.
  * Like：
>>>carbon_rings = Start.choosen_CH(CH_MF) <br>
>>>carbon_rings.to_excel(f"underMCH.xlsx", sheet_name='Sheet1', index=False)
