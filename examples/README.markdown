Factorization in 2025 is one of the great open problems of mathematics and theoretical computer science, marking a turning point in the history of **fundamental science**.

### **Example 1: Factorizing a Single Number**  

#### **Description**  
This example shows how to factorize a single integer using the command line. It is the simplest use case, requiring only a single input number.  

#### **Command**  
```sh
./factor 19984426027
```

#### **Expected Output**  
```
Number: 19984426027
Factors: 2609 (prime), 7659803 (prime)
```

#### **Estimated Execution Time**  
- **Small numbers (< 10 digits)**: Instantaneous (< 1ms)  
- **Medium numbers (10-15 digits)**: A few milliseconds  
- **Larger numbers (15+ digits)**: May take several seconds, depending on system performance  


### **Example 2: Batch Factorization with Formatted Output**  

#### **Description**  
This example demonstrates how to factorize multiple numbers at once using a CSV file for input and outputting the results in JSON format.  

#### **Step 1: Create a CSV File with Numbers**  
Save the following content into a file named `numbers.csv` :  
```csv
number
8724560842712761
1054113196579
16420344313080856219
```

#### **Step 2: Run the Factorization Command**  
```sh
./factor --input-file numbers.csv --output-file results.json --output-json
```

#### **Step 3: Expected Output (JSON Format)**  
The file `results.json` will contain:  
```json
[
    {
        "input": "8724560842712761",
        "factors": [
            {
                "num": "50350709",
                "power": 1,
                "prime": true
            },
            {
                "num": "173275829",
                "power": 1,
                "prime": true
            }
        ]
    },
    {
        "input": "1054113196579",
        "factors": [
            {
                "num": "574061",
                "power": 1,
                "prime": true
            },
            {
                "num": "1836239",
                "power": 1,
                "prime": true
            }
        ]
    },
    {
        "input": "16420344313080856219",
        "factors": [
            {
                "num": "2818298627",
                "power": 1,
                "prime": true
            },
            {
                "num": "5826332297",
                "power": 1,
                "prime": true
            }
        ]
    }
]
```

#### **Estimated Execution Time**  
- **Small numbers (< 15 digits)**: < 10ms per number  
- **Medium numbers (15-30 digits)**: ~50-500ms per number  
- **Larger numbers (30+ digits)**: Seconds to minutes depending on size  

Load the file into a program (Python, R, etc.) to continue in the data pipeline : [![JSON and CSV outputs](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/json-csv-outputs.yaml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/json-csv-outputs.yaml)

### **Example 3: Performance Testing with Verbose Mode**  

#### **Description**  
This example demonstrates how to measure factorization performance using the built-in sample generator and the `--verbose` flag.  

#### **Step 1: Generate Sample Numbers**  
The software includes a built-in number generator for testing. To generate sample numbers ranging from 150-bit to 160-bit :  
```sh
./factor --demand 150 160 --output-file samples.txt --rand-seed 123
```
This will create a file `samples.txt` containing randomly generated numbers based on the seed 123.  

#### **Step 2: Factorize with Verbose Output**  
```sh
./factor --input-file samples.txt --output-csv-extended --verbose 2
```

#### **Example Output**  
```
Input,Factors,Duration (s),Status
1036331975873343534480176146147088186225794603,528033793041353;991506241453181;1979437204147471,0.41,Complete
2630498363192069412551801981938974318484990081,1105375943297653;1248808849566569;1905601384210133,0.34,Complete
5625749556638103914006206263766207107277760581,41486857973238397348399;135603172461675987459019,0.47,Complete
10583771666843650435337231425910227571017879543,18571129;31225147;33182323;60719663;83099501;109009189,0.55,Complete
15205418056784355625162003310385275477360060129,112587861611920715843447;135053795667564377148007,0.55,Complete
25059244536032094373138092233004791182176231637,143712669341;357170154949;677209251977;720899746909,0.70,Complete
46205855979258659499046551995020946278698857633,2843853220513457;4030834075767463;4030834075767463,0.72,Complete
120849154333186312205036763986650288138537680289,276636323879;356553362969;614073425827;1995214398757,0.74,Complete
304917664374965583076514276110184770047563029811,3293978708773679;6584889984856541;14057667705434849,1.68,Complete
546278545875495382441185256230495710955571936601,4421643550154647;8312811386667317;14862179043221299,0.99,Complete
900344327660901048166820504890277547154792935517,6913603345703881;11284236218720999;11540696094952643,1.15,Complete

The results are completely accurate and verified in 8.31 s.
```

#### **Estimated Execution Time**  
Higher `--verbose` levels provide more insights :
- **Verbose level 0 (`--verbose 0`)** → Just factorization, no other messages.  
- **Verbose level 1 (`--verbose 1`)** → Basic progress updates (the default).  
- **Verbose level 2 (`--verbose 2`)** → Status message with execution times.  
- **Verbose level 3 (`--verbose 3`)** → Maintenance messages (Quadratic Sieve).  
- **Verbose level 4 (`--verbose 4`)** → Full debug mode (Quadratic Sieve in-depth analysis).

When factoring numbers larger than 220 bits without a terminal, `--verbose 4` prints one `.` for each percent of progress.

#### **Notes**
- `--demand <bits>` provides a single number of the desired size for testing.  
- The factorization time will vary based on system performance and input size.

### OEIS Data Processing for Factorization

To factorize the OEIS list of choice, Windows users should use **PowerShell**, while macOS and Linux users should use the **Terminal**.

#### Windows (PowerShell Script)
Run `Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope Process` to be able to proceed, then content of `oeis-factoring.ps1` file should be :
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
- Navigate to the directory containing the script using the `cd` command.
- Execute the script using `.\oeis-factoring.ps1` (Windows) or `./oeis-factoring.sh` (macOS and Linux).
- Once downloaded, use the factorization tool with `./factor --input-file oeis-002109.txt --output-csv --verbose 2`.
- This will factorize the entire dataset from OEIS and output the results in CSV format.
- Try OEIS lists, then close the shell (Windows) to ensure any re-application of network rules.

As a workflow : [![Factorize OEIS numbers](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/factorize-oeis.yaml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/factorize-oeis.yaml)

