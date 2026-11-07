# ThermoFinance

Estudiante: Conde Gutierrez Rafael Blas

Aplicación web académica para el desafío **“Aplicación del Método de Gradientes Conjugados Precondicionado en problemas financieros”** de COM-254 Métodos Numéricos.

## Estructura

```text
ThermoFinance/
├── index.html
├── README.md
├── css/
│   └── styles.css
└── js/
    └── app.js
```

## Ejecución local

1. Descarga o clona el repositorio.
2. Abre la carpeta en Visual Studio Code.
3. Ejecuta `index.html` con Live Server, o abre el archivo directamente en un navegador moderno.
4. Se requiere conexión a Internet para cargar Bootstrap, Bootstrap Icons, Google Fonts y Chart.js desde sus CDN.

## Funcionalidad

- Caso de estudio de asignación presupuestaria entre cuatro unidades.
- Tres escenarios: ideal, bajo estrés financiero y mal condicionado.
- Edición de la matriz simétrica A y del vector b.
- Validación de simetría y definición positiva mediante Cholesky.
- Solución numérica de referencia por Cholesky.
- Gradientes Conjugados (GC) y Gradientes Conjugados Precondicionado (GCP) con Jacobi.
- Tolerancia predeterminada `1e-6` y máximo configurable de iteraciones.
- Errores absolutos, relativos y porcentuales por variable.
- Residuo relativo, número de condición de A y condición de la matriz simétrica precondicionada equivalente.
- Gráficas de asignación y convergencia.
- Cuatro pruebas QA exigidas en el documento de desafío.

## Consideraciones matemáticas

La aplicación trabaja con una matriz 4×4, por lo que cumple el requisito mínimo `n >= 4`. La matriz debe ser simétrica y definida positiva. El número de condición en norma 2 se obtiene mediante autovalores de una matriz simétrica, calculados con rotaciones de Jacobi.

Para el sistema precondicionado se calcula el número de condición de `M^(-1/2) A M^(-1/2)`, matriz simétrica similar a `M^(-1) A` cuando `M` es definida positiva. Se utiliza el precondicionador diagonal de Jacobi `M = diag(A)`.

Los escenarios contienen datos sintéticos de demostración. No son datos reales de una empresa ni recomendaciones de inversión.

