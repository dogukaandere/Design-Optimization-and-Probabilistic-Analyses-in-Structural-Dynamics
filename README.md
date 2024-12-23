# Design Optimization and Probabilistic Analyses in Structural Dynamics

This repository contains the MATLAB implementations and documentation for the project **"Design Optimization and Probabilistic Analyses in Structural Dynamics."** The focus is on optimizing a beam profile under specific constraints using deterministic and probabilistic methods. The project tasks are divided into three main areas, each organized in separate folders.

Dieses Repository enthält die MATLAB-Implementierungen und die zugehörige Dokumentation für das Projekt **"Design Optimization and Probabilistic Analyses in Structural Dynamics."** Der Schwerpunkt liegt auf der Optimierung eines Balkenprofils unter Berücksichtigung spezifischer Randbedingungen und der Anwendung deterministischer und probabilistischer Methoden. Die Aufgaben des Projekts sind in drei Hauptbereiche unterteilt, die jeweils in separaten Ordnern organisiert sind.

## Repository Structure / Repository-Struktur

```
Design-Optimization-and-Probabilistic-Analyses-in-Structural-Dynamics/
├── Deterministic-Optimization/
├── Probabilistic-Analysis/
└── Robustness-Optimization/
```

## Task Overview / Aufgabenübersicht

### 1. Deterministic Optimization / Deterministische Optimierung
- **Objective / Ziel**: Minimize the deflection of a beam under specific loading and geometric constraints. / Minimierung der Durchbiegung eines Balkens unter bestimmten Last- und Geometriebedingungen.
- **Approach / Ansatz**: Utilizes MATLAB's `fmincon` function to adjust beam parameters while adhering to boundary conditions. / Verwendung der MATLAB-Funktion `fmincon`, um die Parameter des Balkens anzupassen und dabei die Randbedingungen einzuhalten.
- **Key Highlights / Wichtige Highlights**:
  - Analytical derivative calculations, validated against numerical methods. / Analytische Berechnung der Ableitungen, validiert durch numerische Methoden.
  - Comparison between reference and optimized designs. / Vergleich zwischen Referenz- und optimierten Designs.

### 2. Probabilistic Analysis / Probabilistische Analyse
- **Objective / Ziel**: Evaluate the variability in beam deflection due to stochastic properties of the elasticity modulus. / Bewertung der Variabilität der Balkendurchbiegung aufgrund stochastischer Eigenschaften des Elastizitätsmoduls.
- **Approach / Ansatz**: Applies statistical methods such as Monte Carlo Simulation and First-Order Second Moment (FOSM). / Anwendung statistischer Methoden wie Monte-Carlo-Simulation und First-Order-Second-Moment (FOSM).
- **Key Highlights / Wichtige Highlights**:
  - Distribution fitting and validation using Kolmogorov-Smirnov tests. / Anpassung und Validierung von Wahrscheinlichkeitsverteilungen mit Kolmogorov-Smirnov-Tests.
  - Statistical summaries and visualizations of variability. / Statistische Zusammenfassungen und Visualisierungen der Variabilität.

### 3. Robustness Optimization (Under Revision) / Robustheitsoptimierung (In Überarbeitung)
- **Objective / Ziel**: Simultaneously optimize the mean deflection and its standard deviation for a robust design. / Gleichzeitige Optimierung der mittleren Durchbiegung und ihrer Standardabweichung für ein robustes Design.
- **Approach / Ansatz**: Combines deterministic and probabilistic methods using a multi-objective optimization strategy. / Kombination deterministischer und probabilistischer Methoden unter Verwendung einer Multi-Objective-Optimierungsstrategie.
- **Key Highlights / Wichtige Highlights**:
  - Integration of deterministic insights with probabilistic analyses. / Integration deterministischer Erkenntnisse mit probabilistischen Analysen.
  - Weighting factors to balance deflection minimization and robustness. / Gewichtungsfaktoren zur Balance von Durchbiegungsminimierung und Robustheit.

*Note / Hinweis: Task 3 is currently under revision and will be updated soon. / Aufgabe 3 befindet sich derzeit in Überarbeitung und wird bald aktualisiert.*

---

## Project Highlights / Projekt-Highlights
- **Deterministic Optimization / Deterministische Optimierung**:
  - Adjusts beam geometry to minimize deflection. / Anpassung der Balkengeometrie zur Minimierung der Durchbiegung.
  - Validates analytical derivatives with numerical methods. / Validierung analytischer Ableitungen mit numerischen Methoden.
- **Probabilistic Analysis / Probabilistische Analyse**:
  - Assesses variability in beam deflection using Monte Carlo and FOSM. / Bewertung der Durchbiegungsvariabilität mittels Monte-Carlo-Simulation und FOSM.
  - Validates material property distributions and compares simulation outcomes. / Validierung von Materialeigenschaftsverteilungen und Vergleich der Simulationsergebnisse.
- **Robustness Optimization / Robustheitsoptimierung**:
  - Balances mean deflection and variability for robust designs. / Ausbalancierte Optimierung der mittleren Durchbiegung und ihrer Variabilität.
  - Benchmarks results against deterministic and probabilistic methods. / Benchmarking der Ergebnisse gegen deterministische und probabilistische Methoden.

---

## Navigating the Repository / Navigation im Repository
Each folder contains: / Jeder Ordner enthält:
- MATLAB scripts and functions relevant to the task. / MATLAB-Skripte und Funktionen für die jeweilige Aufgabe.
- A detailed report summarizing methodologies, results, and key insights. / Einen detaillierten Bericht mit Methoden, Ergebnissen und zentralen Erkenntnissen.
- Plots and visualizations illustrating findings. / Plots und Visualisierungen zur Veranschaulichung der Ergebnisse.

### Folder Details / Ordnerdetails
- **Deterministic Optimization / Deterministische Optimierung**:
  - MATLAB code for optimization and derivative validation. / MATLAB-Code für die deterministische Optimierung und die Validierung der Ableitungen.
  - Results of optimized versus reference designs. / Ergebnisse der Referenz- und optimierten Designs.
- **Probabilistic Analysis / Probabilistische Analyse**:
  - Scripts for Monte Carlo simulations and FOSM evaluations. / Skripte für Monte-Carlo-Simulationen und FOSM-Bewertungen.
  - Plots for distribution fitting and validation. / Plots zur Verteilungsanpassung und Validierung.
- **Robustness Optimization / Robustheitsoptimierung**:
  - Scripts for multi-objective optimization with weighting strategies. / Skripte zur Multi-Objective-Optimierung mit Gewichtungsstrategien.
  - Comparative results for deterministic and probabilistic approaches. / Vergleichende Ergebnisse für deterministische und probabilistische Ansätze.

---

## Prerequisites / Voraussetzungen
- **MATLAB Optimization Toolbox** for fmincon. / **MATLAB Optimization Toolbox** für fmincon.

---

## Contact / Kontakt
For further questions or contributions, feel free to contact me via email: / Für weitere Fragen oder Beiträge können Sie mich per E-Mail kontaktieren:
- dogukaan.dere@tuhh.de
- dogu-kaan.dere@outlook.com

---


