## Struktura

```
src/hyppe.py         klient API — samoograniczanie tempa, ponawianie przy 5xx i awariach
                     tunelu Cloudflare
src/judge_scan.py    dyskowy cache porównań; Sędzia jest deterministyczny, więc powtórzone
                     zapytanie nic nie wnosi, a każde trafienie w cache to niezapłacone
                     zapytanie
experiments/         ponumerowane, uruchamiane po kolei; każdy skrypt samodzielny
                     i z punktami kontrolnymi
results/             rankingi, werdykty i finalne sekwencje
```
