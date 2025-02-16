### **Example 1: Factorizing a Single Number (Beginner Level)**  

#### **Description**  
This example shows how to factorize a single integer using the command line. It is the simplest use case, requiring only a single input number.  

#### **Command**  
```sh
factorizer 1234567891011
```

#### **Expected Output**  
```sh
1234567891011 = 2741 × 450283213
```

#### **Estimated Execution Time**  
- **Small numbers (< 10 digits)**: Instantaneous (< 1ms)  
- **Medium numbers (10-15 digits)**: A few milliseconds  
- **Larger numbers (15+ digits)**: May take several seconds, depending on system performance  

#### **Notes**  
- The factorization method used depends on the number’s size.  
- If the number is prime, the output will simply return the number itself.  
- For large numbers (50+ digits), use the **Self-Initializing Quadratic Sieve (SIQS)** method for better efficiency.  


Super ! Pour le niveau intermédiaire, on va inclure :  

- **Batch processing** : Factorisation de plusieurs nombres en une seule commande.  
- **Fichiers CSV & JSON** : Entrées et sorties formatées.  
- **Exemples de commandes avec estimation du temps d’exécution**.  


### **Example 2: Batch Factorization with Formatted Output (Intermediate Level)**  

#### **Description**  
This example demonstrates how to factorize multiple numbers at once using a CSV file for input and outputting the results in JSON format.  

#### **Step 1: Create a CSV File with Numbers**  
Save the following content into a file named `numbers.csv`:  
```csv
number
1234567891011
987654321001
112233445566778899
```

#### **Step 2: Run the Factorization Command**  
```sh
factorizer --input numbers.csv --output results.json
```

#### **Step 3: Expected Output (JSON Format)**  
The file `results.json` will contain:  
```json
[
  { "number": 1234567891011, "factors": [2741, 450283213] },
  { "number": 987654321001, "factors": [17, 57945136529] },
  { "number": 112233445566778899, "factors": [19, 59, 97, 104729, 16769023] }
]
```

#### **Estimated Execution Time**  
- **Small numbers (< 15 digits)**: < 10ms per number  
- **Medium numbers (15-30 digits)**: ~50-500ms per number  
- **Larger numbers (30+ digits)**: Seconds to minutes depending on size  

#### **Notes**  
- Supports CSV and JSON formats for both input and output.  
- The software automatically selects the best factorization method.  
- For large numbers (50+ digits), consider enabling SIQS for better performance.

### **Example 3: Performance Testing with Verbose Mode (Advanced Users)**  

#### **Description**  
This example demonstrates how to measure factorization performance using the built-in sample generator and the `--verbose` flag.  

#### **Step 1: Generate Sample Numbers**  
The software includes a built-in number generator for testing. To generate 5 sample numbers:  
```sh
factorizer --generate 5 --output samples.txt
```
This will create a file `samples.txt` containing randomly generated numbers of various sizes.  

#### **Step 2: Factorize with Verbose Output**  
```sh
factorizer --input samples.txt --verbose 2
```

#### **Example Output**  
```sh
[INFO] Factoring 123456789101112131415...
[INFO] Using Quadratic Sieve
[INFO] Step 1: Sieving completed in 1.23s
[INFO] Step 2: Matrix reduction completed in 0.45s
[INFO] Step 3: Root extraction completed in 0.02s
123456789101112131415 = 3571 × 345678902341

[INFO] Factoring 987654321987654321...
[INFO] Using Pollard’s Rho (small number detected)
[INFO] Factorization completed in 0.0005s
987654321987654321 = 3 × 17 × 193 × 17719 × 269741
```

#### **Estimated Execution Time**  
- **Verbose level 1 (`--verbose 1`)** → Basic progress updates (~minimal impact on speed).  
- **Verbose level 2 (`--verbose 2`)** → Detailed step-by-step execution times (~2-5% overhead).  
- **Verbose level 3 (`--verbose 3`)** → Full debug mode (~significant slowdown, useful for in-depth analysis).  

#### **Notes**  
- `--generate <N>` provides random numbers for testing.  
- Higher `--verbose` levels provide more insights but slow down execution.  
- The factorization time will vary based on system performance and input size.  

### OEIS Data Processing for Factorization

 Windows users should use **PowerShell**, while macOS and Linux users should use the **Terminal**.

#### Windows (PowerShell Script)
Run `Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser` to be able to process. The contents of the file `oeis-factoring.ps1` is :
```powershell
$oeis = "002109"  # Select the Hyperfactorials from OEIS
$url = "https://oeis.org/A$oeis/b$oeis.txt"
$output = "oeis-$oeis.txt"

# Download the file
Invoke-WebRequest -Uri $url -OutFile "oeis-raw.txt"

# Prune the file
Get-Content "oeis-raw.txt" | ForEach-Object {
    if ($_ -match "^\s*\d+\s+(\S+)") { $matches[1] }
} | Set-Content $output

Write-Host "The file is ready for factoring at $output"
```

#### macOS (Bash Script)
Run `chmod 0700 oeis-factoring.sh` after creating the file `oeis-factoring.sh` with :
```bash
#!/bin/bash

oeis="002109"  # Select the Hyperfactorials from OEIS
url="https://oeis.org/A$oeis/b$oeis.txt"
output="oeis-$oeis.txt"

# Download the file
curl -s "$url" -o "oeis-raw.txt"

# Prune the file
awk '{if ($1 ~ /^[0-9]+$/) print $2}' oeis-raw.txt > "$output"

echo "The file is ready for factoring at $output"
```

#### Linux (Bash Script)
Run `chmod 0700 oeis-factoring.sh` after creating the file `oeis-factoring.sh` with :
```bash
#!/bin/bash

oeis="002109"  # Select the Hyperfactorials from OEIS
url="https://oeis.org/A$oeis/b$oeis.txt"
output="oeis-$oeis.txt"

# Download the file
wget -q "$url" -O "oeis-raw.txt"

# Prune the file
awk '{if ($1 ~ /^[0-9]+$/) print $2}' oeis-raw.txt > "$output"

echo "The file is ready for factoring at $output"
```

#### Usage
Once downloaded, use the factorization tool with `./factor --input-file oeis-002109.txt --output-csv --verbose 2`. This will factorize the entire dataset from OEIS and output the results in CSV format.
