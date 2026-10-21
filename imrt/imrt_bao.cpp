#include <functional>
#include "imrt_bao.h"

#include "../emilibase.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace emili {
namespace imrt {

//=== Representación de solución ===============================================

emili::Solution* BaoSolution::clone()
{
    // Crea un nuevo objeto en heap con los mismos ángulos activos e intensidades
    BaoSolution* s = new BaoSolution(active_angles_, (int)intensities_.size());

    // Copia los vectores por valor — el clon es independiente del original
    s->angle_degrees_    = angle_degrees_;
    s->intensities_      = intensities_;
    s->ptv_underdose_sq_ = ptv_underdose_sq_;
    s->oar_overdose_sq_  = oar_overdose_sq_;

    // Copia el valor objetivo para no pagar otra evaluación FMO
    s->solution_value = solution_value;

    return s;
}

std::string BaoSolution::getSolutionRepresentation()
{
    std::ostringstream oss;

    // Si no fue evaluada aún, muestra los índices activos (angle_degrees_ vacío)
    const std::vector<int>& deg = angle_degrees_.empty() ? active_angles_ : angle_degrees_;

    oss << "angles=[";
    for (int i = 0; i < (int)deg.size(); ++i) {
        if (i) oss << ",";   // separador de coma entre elementos (omite antes del primero)
        oss << deg[i];
    }
    oss << "]" << (angle_degrees_.empty() ? "(idx)" : "deg");

    // Agrega el valor objetivo con precisión fija de 2 decimales
    oss << " f=" << std::fixed << std::setprecision(2) << solution_value;

    return oss.str();
}


//=== Evaluación y registro ====================================================

void BaoProblem::openCsvLog(const std::string& path)
{
    csv_file_.open(path);
    if (csv_file_.is_open()) {
        // `cached` distingue una visita de una resolución FMO nueva. Después
        // vienen dos columnas por órgano: ptv_<nombre>_u2 (subdosis al
        // cuadrado, el término que aporta ese PTV a w_under*sum(u^2)) y
        // oar_<nombre>_v2 (sobredosis al cuadrado, término de w_over*sum(v^2))
        // — el desglose que el objetivo agregado no muestra por sí solo.
        csv_file_ << "eval,angles_deg,objective,cached";
        for (const auto& name : ptvOrganNames()) csv_file_ << ",ptv_" << name << "_u2";
        for (const auto& name : oarOrganNames()) csv_file_ << ",oar_" << name << "_v2";
        csv_file_ << "\n";
        std::cout << "  CSV log: " << path << "\n";
    } else {
        // Falla silenciosa: el algoritmo continúa sin log si el archivo no se puede abrir
        std::cerr << "[BAO] Could not open CSV log: " << path << "\n";
    }
}

double BaoProblem::calcObjectiveFunctionValue(emili::Solution& s)
{
    // Downcast seguro: s siempre es BaoSolution en este contexto
    BaoSolution& bs = static_cast<BaoSolution&>(s);

    // El FMO depende del conjunto de ángulos, no del orden de los slots.
    // Canonicalizar la clave permite reutilizar también permutaciones equivalentes.
    std::vector<int> cache_key = bs.active_angles_;
    std::sort(cache_key.begin(), cache_key.end());

    auto cached = fmo_cache_.find(cache_key);
    if (cached != fmo_cache_.end()) {
        bs.intensities_      = cached->second.intensities;
        bs.ptv_underdose_sq_ = cached->second.ptv_underdose_sq;
        bs.oar_overdose_sq_  = cached->second.oar_overdose_sq;
        last_eval_cached_ = true;
        return cached->second.objective;
    }

    // Llamada más costosa del sistema: resuelve el FMO con los ángulos activos.
    FmoResult res = fmo_.solve(bs.active_angles_);

    // Guarda objetivo, intensidades y desglose por órgano para evitar
    // futuras llamadas idénticas al solver FMO.
    CachedFmoResult entry;
    entry.intensities      = std::move(res.intensities);
    entry.objective        = res.objective;
    entry.ptv_underdose_sq = std::move(res.ptv_underdose_sq);
    entry.oar_overdose_sq  = std::move(res.oar_overdose_sq);
    auto inserted = fmo_cache_.insert(std::make_pair(cache_key, std::move(entry)));

    bs.intensities_      = inserted.first->second.intensities;
    bs.ptv_underdose_sq_ = inserted.first->second.ptv_underdose_sq;
    bs.oar_overdose_sq_  = inserted.first->second.oar_overdose_sq;
    last_eval_cached_ = false;

    // Devuelve solo el valor numérico del objetivo; las intensidades y el
    // desglose por órgano quedan en bs
    return inserted.first->second.objective;
}

double BaoProblem::evaluateSolution(emili::Solution& s)
{
    // Resuelve el FMO y llena bs.intensities_ como efecto secundario
    double f = calcObjectiveFunctionValue(s);

    // Persiste el valor en la clase base para que el algoritmo pueda comparar soluciones
    s.setSolutionValue(f);

    // Downcast seguro: necesario para acceder a angle_degrees_ y active_angles_
    BaoSolution& bs = static_cast<BaoSolution&>(s);

    // Convierte índices activos a grados reales del catálogo — llena angle_degrees_
    bs.angle_degrees_.resize(bs.active_angles_.size());
    for (int i = 0; i < (int)bs.active_angles_.size(); ++i)
        bs.angle_degrees_[i] = angleDegree(bs.active_angles_[i]);

    // Logging en consola solo si verbose_ está activo — no afecta el resultado
    if (verbose_) {
        std::cout << "  BAO eval: angles=[";
        for (int i = 0; i < (int)bs.angle_degrees_.size(); ++i) {
            if (i) std::cout << ",";
            std::cout << bs.angle_degrees_[i];
        }
        std::cout << "deg] -> f=" << std::fixed << std::setprecision(2) << f;
        if (last_eval_cached_) std::cout << " [cached]";
        std::cout << "\n";
    }

    if (csv_file_.is_open()) {
        // Preincremento: el contador arranca en 1
        csv_file_ << ++eval_count_ << ",\"";
        for (int i = 0; i < (int)bs.angle_degrees_.size(); ++i) {
            if (i) csv_file_ << ";";   // separador interno ; para no romper el formato CSV
            csv_file_ << bs.angle_degrees_[i];
        }
        // 6 decimales en CSV — más precisión que el display en consola
        csv_file_ << "\"," << std::fixed << std::setprecision(6) << f
                  << "," << (last_eval_cached_ ? "true" : "false");
        // Desglose por órgano, mismo orden que el encabezado armó en openCsvLog
        for (double u2 : bs.ptv_underdose_sq_) csv_file_ << "," << std::setprecision(6) << u2;
        for (double v2 : bs.oar_overdose_sq_)  csv_file_ << "," << std::setprecision(6) << v2;
        csv_file_ << "\n";
        // flush fuerza escritura al disco — útil si el proceso se interrumpe
        csv_file_.flush();
    }

    return f;
}

//=== Soluciones iniciales =====================================================

emili::Solution* FirstKAnglesInit::generateEmptySolution()
{
    // Ordena índices, no el catálogo: bao_.angleDegree() mantiene su
    // correspondencia con las matrices de dosis cargadas para cada ángulo.
    // Filtra a múltiplos de 5 por el mismo motivo que RandomKAnglesInit: si
    // el punto de partida no cae en la reja de 5° que usa el vecindario de
    // shift, el espacio alcanzable queda desfasado según la estrategia de
    // inicialización usada, y las comparaciones dejan de ser sobre el mismo
    // espacio de búsqueda. Si restrict_to_5deg_grid_ es false (catálogo
    // "unrestricted"), no se filtra: se usan todos los índices del catálogo.
    std::vector<int> degree_order;
    for (int i = 0; i < bao_.nAngles(); ++i) {
        if (!restrict_to_5deg_grid_ || bao_.angleDegree(i) % 5 == 0) degree_order.push_back(i);
    }
    std::sort(degree_order.begin(), degree_order.end(),
        [this](int a, int b) { return bao_.angleDegree(a) < bao_.angleDegree(b); });

    // Selecciona los K ángulos de menor grado real (ya filtrados a múltiplos
    // de 5). Por ejemplo, en el catálogo de CERR_Prostate produce
    // 0°, 5°, 10°, 15° para K=4.
    std::vector<int> angles(degree_order.begin(),
                            degree_order.begin() + bao_.K());
    std::sort(angles.begin(), angles.end());

    // Crea la solución en heap — sin evaluar: angle_degrees_ vacío, solution_value = 0
    return new BaoSolution(angles, bao_.nDimlets());
}

emili::Solution* FirstKAnglesInit::generateSolution()
{
    // Construye la estructura vacía con los primeros K ángulos
    emili::Solution* s = generateEmptySolution();

    // Evalúa: resuelve el FMO y llena angle_degrees_ y solution_value
    bao_.evaluateSolution(*s);
    return s;
}

emili::Solution* RandomKAnglesInit::generateEmptySolution()
{
    int K = bao_.K();  // cuántos ángulos debe tener la solución

    // Restringe el pool a índices cuyo grado real es múltiplo de 5: si el
    // vecindario de shift avanza en pasos múltiplos de 5, la solución inicial
    // debe partir en esa misma reja, o el espacio alcanzable por shift queda
    // desfasado según el punto de partida (y las comparaciones entre corridas
    // dejan de ser sobre el mismo espacio de búsqueda). Si
    // restrict_to_5deg_grid_ es false (catálogo "unrestricted"), no se
    // filtra: el pool son todos los índices del catálogo.
    std::vector<int> pool;
    for (int i = 0; i < bao_.nAngles(); ++i) {
        if (!restrict_to_5deg_grid_ || bao_.angleDegree(i) % 5 == 0) pool.push_back(i);
    }
    int n = (int)pool.size();

    // Fisher-Yates parcial sobre el pool filtrado: mezcla solo los primeros K
    // elementos — O(K) en vez de O(n). No se usa std::shuffle porque
    // emili::generateRandomNumber() no es compatible con std.
    for (int i = 0; i < K; ++i) {
        int j = i + emili::generateRandomNumber() % (n - i);  // índice aleatorio en [i, n-1]
        std::swap(pool[i], pool[j]);
    }

    // Toma los primeros K elementos mezclados como ángulos activos
    std::vector<int> angles(pool.begin(), pool.begin() + K);

    // Ordena los índices para mantener consistencia con el resto del algoritmo
    std::sort(angles.begin(), angles.end());

    // Crea la solución en heap — sin evaluar: angle_degrees_ vacío, solution_value = 0
    return new BaoSolution(angles, bao_.nDimlets());
}

emili::Solution* RandomKAnglesInit::generateSolution()
{
    // Construye la estructura vacía con K ángulos aleatorios
    emili::Solution* s = generateEmptySolution();

    // Evalúa: resuelve el FMO y llena angle_degrees_ y solution_value
    bao_.evaluateSolution(*s);
    return s;
}

//=== Vecindario: swap de ángulos ==============================================

int AngleSwapNeighborhood::size()
{
    int K  = bao_.K();                      // ángulos activos
    int na = bao_.nAngles();                // total de ángulos en el catálogo

    // Vecindario = K opciones para sacar × (na - K) opciones para meter
    // Ej: K=5, na=36 → 5 × 31 = 155 vecinos
    return K * (na - K);
}

emili::Neighborhood::NeighborhoodIterator
AngleSwapNeighborhood::begin(emili::Solution* base)
{
    // Downcast seguro: necesitamos active_angles_ que no existe en Solution base
    BaoSolution* bs = static_cast<BaoSolution*>(base);

    // Guarda copia del estado base — computeStep siempre parte de acá para cada swap
    base_angles_ = bs->active_angles_;

    // Marca qué ángulos están activos en un vector de flags
    std::vector<bool> active_flag(n_angles_, false);
    for (int ai : base_angles_) active_flag[ai] = true;

    // Recolecta todos los ángulos que NO están activos — candidatos para entrar en un swap
    inactive_list_.clear();
    for (int a = 0; a < n_angles_; ++a)
        if (!active_flag[a]) inactive_list_.push_back(a);

    // Inicializa cursores al inicio del recorrido
    cur_active_idx_   = 0;    // primer ángulo activo a sacar
    cur_inactive_idx_ = 0;    // primer ángulo inactivo a meter
    first_            = true; // evita avanzar cursores en la primera llamada a computeStep

    // Devuelve el iterador listo para que el framework llame a computeStep K×(n-K) veces
    return emili::Neighborhood::NeighborhoodIterator(this, base);
}

void AngleSwapNeighborhood::reset()
{
    // Reinicia el recorrido al primer vecino — usado cuando el framework relanza la búsqueda
    cur_active_idx_   = 0;
    cur_inactive_idx_ = 0;
    first_            = true;
}

emili::Solution* AngleSwapNeighborhood::computeStep(emili::Solution* step)
{
    if (!first_) {
        // Avanza al siguiente inactivo; si se agotaron, pasa al siguiente activo
        ++cur_inactive_idx_;
        if (cur_inactive_idx_ >= (int)inactive_list_.size()) {
            ++cur_active_idx_;
            cur_inactive_idx_ = 0;
        }
    }
    first_ = false;

    // Fin del vecindario: se recorrieron todos los swaps posibles
    if (cur_active_idx_ >= (int)base_angles_.size()) return nullptr;
    if (inactive_list_.empty())                       return nullptr;

    // Aplica el swap: reemplaza el ángulo activo actual por el inactivo actual
    BaoSolution* bs = static_cast<BaoSolution*>(step);
    bs->active_angles_ = base_angles_;
    bs->active_angles_[cur_active_idx_] = inactive_list_[cur_inactive_idx_];
    std::sort(bs->active_angles_.begin(), bs->active_angles_.end());

    // Evalúa el vecino generado
    bao_.evaluateSolution(*bs);
    return bs;
}

void AngleSwapNeighborhood::reverseLastMove(emili::Solution* step)
{
    // Restaura los ángulos activos al estado base — el framework ya resetea solution_value
    // Las intensidades quedan desactualizadas pero se sobreescriben en el próximo computeStep
    BaoSolution* bs = static_cast<BaoSolution*>(step);
    bs->active_angles_ = base_angles_;
}

emili::Solution* AngleSwapNeighborhood::random(emili::Solution* s)
{
    BaoSolution* bs = static_cast<BaoSolution*>(s);

    // Si no hay inactivos no se puede hacer swap — devuelve clon sin cambios
    if (inactive_list_.empty()) return s->clone();

    // Elige un ángulo activo e inactivo al azar
    int ai = emili::generateRandomNumber() % bs->active_angles_.size();
    int ii = emili::generateRandomNumber() % inactive_list_.size();

    // Aplica el swap sobre una copia para no mutar la solución original
    BaoSolution* nb = static_cast<BaoSolution*>(s->clone());
    nb->active_angles_[ai] = inactive_list_[ii];
    std::sort(nb->active_angles_.begin(), nb->active_angles_.end());
    bao_.evaluateSolution(*nb);
    return nb;
}

//=== Vecindario: shift de ángulos ==============================================

void AngleShiftNeighborhood::buildDegreeOrder()
{
    // Construye la permutación de índices de catálogo ordenados por grado real ascendente
    // (el array crudo del catálogo puede estar en orden lexicográfico de string, no numérico)
    degree_order_.resize(n_angles_);
    for (int i = 0; i < n_angles_; ++i) degree_order_[i] = i;
    std::sort(degree_order_.begin(), degree_order_.end(),
        [this](int a, int b) { return bao_.angleDegree(a) < bao_.angleDegree(b); });

    // Lookup inverso: índice de catálogo -> su posición en degree_order_
    degree_rank_.resize(n_angles_);
    for (int pos = 0; pos < n_angles_; ++pos) degree_rank_[degree_order_[pos]] = pos;
}

int AngleShiftNeighborhood::size()
{
    // Cota superior: 2 movimientos (±step) por cada ángulo activo.
    // El conteo real puede ser menor por colisiones con otros ángulos activos.
    return 2 * bao_.K();
}

emili::Neighborhood::NeighborhoodIterator
AngleShiftNeighborhood::begin(emili::Solution* base)
{
    // First/Best Improvement llama a begin() al iniciar cada ronda de búsqueda.
    // Guardamos la solución actual porque TODOS los vecinos de esta ronda deben
    // construirse a partir de la misma base.
    BaoSolution* bs = static_cast<BaoSolution*>(base);

    // Modo circular: antes de sobreescribir base_angles_ (que todavía guarda
    // los ángulos de la RONDA ANTERIOR en este punto) con la nueva base,
    // comparamos ambos vectores slot a slot para encontrar cuál slot cambió
    // por el último movimiento aceptado. Ese slot es el punto de arranque de
    // la rotación de esta ronda. Si no hay base_angles_ previo (primera
    // llamada) o random_order_ tiene prioridad, arranca en el slot 0 (igual
    // que los otros modos).
    int rotation_start = 0;
    if (circular_order_ && !random_order_ && !base_angles_.empty()) {
        const std::vector<int>& new_angles = bs->active_angles_;
        for (size_t i = 0; i < base_angles_.size() && i < new_angles.size(); ++i) {
            if (base_angles_[i] != new_angles[i]) { rotation_start = (int)i; break; }
        }
    }

    base_angles_ = bs->active_angles_;

    // El catálogo puede estar almacenado como 0,100,10,110,...; este lookup
    // permite que una posición represente el siguiente grado real: 0,10,20,...
    buildDegreeOrder();

    // Orden de generación (por defecto, random_order_ = false y
    // circular_order_ = false):
    //   slot 0: -step, +step
    //   slot 1: -step, +step
    //   ...
    // Con random_order_ = true, slot_order_ es una permutación aleatoria de
    // 0..K-1 (Fisher-Yates, mismo estilo que AngleShiftMultiPerturbation::
    // perturb) y se recorre esa en vez del orden fijo. Con circular_order_ =
    // true (y random_order_ = false), slot_order_ es una rotación de 0..K-1
    // que arranca en rotation_start (el slot recién modificado por la ronda
    // anterior) en vez de reiniciar siempre en el slot 0 — la teoría es que,
    // en un paisaje continuo (ángulos irradiando un volumen fijo), seguir
    // explorando alrededor del mismo ángulo que acaba de mejorar es más
    // prometedor que volver siempre al slot 0. En ambos casos solo cambia el
    // ORDEN en que First Improvement descubre vecinos; Best Improvement
    // evalúa el vecindario completo igual, así que no le afecta salvo en
    // desempates.
    int K = (int)base_angles_.size();
    slot_order_.resize(K);
    if (random_order_) {
        for (int i = 0; i < K; ++i) slot_order_[i] = i;
        for (int i = K - 1; i > 0; --i) {
            int j = emili::generateRandomNumber() % (i + 1);
            std::swap(slot_order_[i], slot_order_[j]);
        }
    } else if (circular_order_) {
        for (int i = 0; i < K; ++i) slot_order_[i] = (rotation_start + i) % K;
    } else {
        for (int i = 0; i < K; ++i) slot_order_[i] = i;
    }

    cur_active_idx_ = 0;
    cur_dir_        = 0;    // 0: resta; 1: suma
    first_          = true;

    return emili::Neighborhood::NeighborhoodIterator(this, base);
}

void AngleShiftNeighborhood::reset()
{
    // Reinicia el recorrido al primer vecino — usado cuando el framework relanza la búsqueda.
    // slot_order_ (fijo o aleatorio) NO se regenera aquí: se conserva de begin(),
    // igual que base_angles_.
    cur_active_idx_ = 0;
    cur_dir_        = 0;
    first_          = true;
}

emili::Solution* AngleShiftNeighborhood::computeStep(emili::Solution* step)
{
    // Esta es la función usada por First Improvement y Best Improvement.
    // Genera UN vecino determinista por llamada — la única aleatoriedad de
    // esta clase entra en begin() (sorteo de slot_order_ si random_order_),
    // no aquí: dado un slot_order_ ya fijado, el recorrido es determinista.
    BaoSolution* bs = static_cast<BaoSolution*>(step);

    // step_ cuenta posiciones del catálogo ordenado por grados. En una instancia
    // con ángulos cada 10°, step_=1 equivale a desplazar exactamente 10°.
    // La normalización también permite desplazamientos circulares.
    int step_mod = ((step_ % n_angles_) + n_angles_) % n_angles_;

    while (true) {
        if (!first_) {
            // Después de probar -step, prueba +step sobre EL MISMO slot.
            // Cuando ambas direcciones terminan, avanza al siguiente slot.
            ++cur_dir_;
            if (cur_dir_ >= 2) {
                ++cur_active_idx_;
                cur_dir_ = 0;
            }
        }
        first_ = false;

        // No quedan combinaciones (slot, dirección): termina el iterador.
        if (cur_active_idx_ >= (int)base_angles_.size()) return nullptr;

        // slot_order_[cur_active_idx_] es el slot real a desplazar: identidad
        // en modo fijo (slot == cur_active_idx_), o la permutación sorteada en
        // begin() en modo random_order_.
        int slot = slot_order_[cur_active_idx_];

        // Obtiene el ángulo que ocupa el slot actual y busca su posición dentro
        // del catálogo ordenado numéricamente.
        int active_catalog_idx = base_angles_[slot];
        int pos     = degree_rank_[active_catalog_idx];

        // Aplica -step o +step con wrap-around circular. Ejemplos para step_=1:
        //   20° - 10° = 10°
        //    0° - 10° = 350°
        //  350° + 10° = 0°
        int new_pos = (cur_dir_ == 0)
            ? (pos - step_mod + n_angles_) % n_angles_
            : (pos + step_mod) % n_angles_;
        int candidate = degree_order_[new_pos];

        // Una solución BAO no puede contener dos veces el mismo ángulo. Si el
        // candidato ya ocupa otro slot, esta combinación se omite y el while
        // continúa automáticamente con la dirección o el slot siguiente.
        bool collision = false;
        for (int a : base_angles_) {
            if (a == candidate) { collision = true; break; }
        }
        if (collision) continue;

        // Siempre reconstruye el vecino desde base_angles_. Solo reemplaza el
        // slot actual y NO ordena el vector, preservando la identidad del slot:
        //   base   [20,120,230,340]
        //   vecino [10,120,230,340]
        bs->active_angles_ = base_angles_;
        bs->active_angles_[slot] = candidate;

        // evaluateSolution consulta primero la caché FMO. El solver FMO solo se
        // ejecuta si este conjunto de ángulos todavía no fue evaluado.
        bao_.evaluateSolution(*bs);
        return bs;
    }
}

void AngleShiftNeighborhood::reverseLastMove(emili::Solution* step)
{
    // Restaura los ángulos activos al estado base — el framework ya resetea solution_value
    BaoSolution* bs = static_cast<BaoSolution*>(step);
    bs->active_angles_ = base_angles_;
}

emili::Solution* AngleShiftNeighborhood::random(emili::Solution* s)
{
    // Esta función NO participa en `first ... nangshift` ni en
    // `best ... nangshift`. Solo se usa si otro componente pide explícitamente
    // un vecino aleatorio al Neighborhood (por ejemplo, una perturbación).
    BaoSolution* bs = static_cast<BaoSolution*>(s);

    // Si random() se llama sin haber pasado antes por begin(), construye el orden por grado
    if ((int)degree_order_.size() != n_angles_) buildDegreeOrder();

    // Marca qué ángulos están activos, para detectar colisiones
    std::vector<bool> active_flag(n_angles_, false);
    for (int a : bs->active_angles_) active_flag[a] = true;

    // step_ define la escala de la perturbación, no su valor exacto: un
    // desplazamiento de magnitud constante preserva el residuo módulo step_
    // del ángulo durante TODA la corrida de ILS si coincide con el step del
    // vecindario de búsqueda local (nangshift) — la perturbación quedaría
    // atrapada en la misma clase módulo step_ que le tocó al ángulo en la
    // inicialización, sin poder escapar nunca de ella. Además, una magnitud
    // menor o igual a step_ es un desplazamiento débil: cae dentro (o al
    // borde) del alcance de un solo paso de la búsqueda local, que puede
    // reconverger casi de inmediato al mismo óptimo. Por eso la magnitud se
    // sortea en (step_, 2*step_): siempre estrictamente mayor que un paso de
    // búsqueda local (sale genuinamente de ese vecindario) y nunca múltiplo
    // de step_ (rompe la invariante de residuo, garantizado).
    int base_step = ((step_ % n_angles_) + n_angles_) % n_angles_;
    if (base_step < 1) base_step = 1;
    int span = base_step - 1;   // cuántos valores hay en (base_step, 2*base_step)
    if (span < 1) span = 1;     // step_ muy chico (1): no hay hueco real, degrada a [1, n_angles_-1]

    // Reintento acotado: prueba unas pocas combinaciones (ángulo, dirección, magnitud) al azar
    const int max_attempts = 8;
    for (int attempt = 0; attempt < max_attempts; ++attempt) {
        int ai       = emili::generateRandomNumber() % bs->active_angles_.size();
        int dir      = emili::generateRandomNumber() % 2;
        int step_mod = (base_step > 1)
            ? base_step + 1 + (emili::generateRandomNumber() % span)
            : 1 + (emili::generateRandomNumber() % (n_angles_ - 1));

        int active_catalog_idx = bs->active_angles_[ai];
        int pos     = degree_rank_[active_catalog_idx];
        int new_pos = (dir == 0)
            ? (pos - step_mod + n_angles_) % n_angles_
            : (pos + step_mod) % n_angles_;
        int candidate = degree_order_[new_pos];

        if (active_flag[candidate]) continue;   // colisión: reintenta con otra combinación

        // Aplica el shift sobre una copia para no mutar la solución original
        // y conserva el slot elegido para que la trayectoria identifique
        // inequívocamente qué ángulo fue desplazado.
        BaoSolution* nb = static_cast<BaoSolution*>(s->clone());
        nb->active_angles_[ai] = candidate;
        bao_.evaluateSolution(*nb);
        return nb;
    }

    // No se encontró un movimiento válido tras los reintentos: devuelve clon sin cambios
    return s->clone();
}

//=== Vecindario: shift de ángulos, múltiples steps (versión inclusiva) =======

void AngleMultiShiftNeighborhood::buildDegreeOrder()
{
    // Idéntico a AngleShiftNeighborhood::buildDegreeOrder — ver ese comentario.
    degree_order_.resize(n_angles_);
    for (int i = 0; i < n_angles_; ++i) degree_order_[i] = i;
    std::sort(degree_order_.begin(), degree_order_.end(),
        [this](int a, int b) { return bao_.angleDegree(a) < bao_.angleDegree(b); });

    degree_rank_.resize(n_angles_);
    for (int pos = 0; pos < n_angles_; ++pos) degree_rank_[degree_order_[pos]] = pos;
}

int AngleMultiShiftNeighborhood::size()
{
    // Cota superior: 2 movimientos (±step) por cada ángulo activo, por cada step
    // de la lista. El conteo real puede ser menor por colisiones.
    return 2 * bao_.K() * (int)steps_.size();
}

emili::Neighborhood::NeighborhoodIterator
AngleMultiShiftNeighborhood::begin(emili::Solution* base)
{
    BaoSolution* bs = static_cast<BaoSolution*>(base);
    base_angles_ = bs->active_angles_;

    buildDegreeOrder();

    // Orden de generación:
    //   slot 0, step steps_[0]: -step, +step
    //   slot 0, step steps_[1]: -step, +step
    //   ...
    //   slot 1, step steps_[0]: -step, +step
    //   ...
    cur_active_idx_ = 0;
    cur_step_idx_   = 0;
    cur_dir_        = 0;    // 0: resta; 1: suma
    first_          = true;

    return emili::Neighborhood::NeighborhoodIterator(this, base);
}

void AngleMultiShiftNeighborhood::reset()
{
    cur_active_idx_ = 0;
    cur_step_idx_   = 0;
    cur_dir_        = 0;
    first_          = true;
}

emili::Solution* AngleMultiShiftNeighborhood::computeStep(emili::Solution* step)
{
    // Mismo esqueleto que AngleShiftNeighborhood::computeStep, con una
    // dimensión extra de avance: dirección -> step de la lista -> slot activo.
    BaoSolution* bs = static_cast<BaoSolution*>(step);

    while (true) {
        if (!first_) {
            ++cur_dir_;
            if (cur_dir_ >= 2) {
                cur_dir_ = 0;
                ++cur_step_idx_;
                if (cur_step_idx_ >= (int)steps_.size()) {
                    cur_step_idx_ = 0;
                    ++cur_active_idx_;
                }
            }
        }
        first_ = false;

        // No quedan combinaciones (slot, step, dirección): termina el iterador.
        if (cur_active_idx_ >= (int)base_angles_.size()) return nullptr;

        int active_catalog_idx = base_angles_[cur_active_idx_];
        int pos = degree_rank_[active_catalog_idx];

        int step_mod = ((steps_[cur_step_idx_] % n_angles_) + n_angles_) % n_angles_;
        int new_pos = (cur_dir_ == 0)
            ? (pos - step_mod + n_angles_) % n_angles_
            : (pos + step_mod) % n_angles_;
        int candidate = degree_order_[new_pos];

        // Colisión con otro slot activo, o con el propio ángulo del slot
        // (puede pasar si dos steps distintos generan el mismo candidato):
        // se omite y el while sigue con la siguiente combinación.
        bool collision = false;
        for (int a : base_angles_) {
            if (a == candidate) { collision = true; break; }
        }
        if (collision) continue;

        bs->active_angles_ = base_angles_;
        bs->active_angles_[cur_active_idx_] = candidate;

        bao_.evaluateSolution(*bs);
        return bs;
    }
}

void AngleMultiShiftNeighborhood::reverseLastMove(emili::Solution* step)
{
    BaoSolution* bs = static_cast<BaoSolution*>(step);
    bs->active_angles_ = base_angles_;
}

emili::Solution* AngleMultiShiftNeighborhood::random(emili::Solution* s)
{
    // No participa en `first ... nangshiftmulti` ni en `best ... nangshiftmulti`
    // (igual que AngleShiftNeighborhood::random) — ver ese comentario.
    BaoSolution* bs = static_cast<BaoSolution*>(s);

    if ((int)degree_order_.size() != n_angles_) buildDegreeOrder();

    std::vector<bool> active_flag(n_angles_, false);
    for (int a : bs->active_angles_) active_flag[a] = true;

    const int max_attempts = 8;
    for (int attempt = 0; attempt < max_attempts; ++attempt) {
        int ai       = emili::generateRandomNumber() % bs->active_angles_.size();
        int step_idx = emili::generateRandomNumber() % steps_.size();
        int dir      = emili::generateRandomNumber() % 2;

        int active_catalog_idx = bs->active_angles_[ai];
        int pos       = degree_rank_[active_catalog_idx];
        int step_mod  = ((steps_[step_idx] % n_angles_) + n_angles_) % n_angles_;
        int new_pos   = (dir == 0)
            ? (pos - step_mod + n_angles_) % n_angles_
            : (pos + step_mod) % n_angles_;
        int candidate = degree_order_[new_pos];

        if (active_flag[candidate]) continue;

        BaoSolution* nb = static_cast<BaoSolution*>(s->clone());
        nb->active_angles_[ai] = candidate;
        bao_.evaluateSolution(*nb);
        return nb;
    }

    return s->clone();
}

//=== Perturbación: shift multi-ángulo, slots distintos ========================

void AngleShiftMultiPerturbation::buildDegreeOrder()
{
    // Idéntico a AngleShiftNeighborhood::buildDegreeOrder — ver ese comentario.
    degree_order_.resize(n_angles_);
    for (int i = 0; i < n_angles_; ++i) degree_order_[i] = i;
    std::sort(degree_order_.begin(), degree_order_.end(),
        [this](int a, int b) { return bao_.angleDegree(a) < bao_.angleDegree(b); });

    degree_rank_.resize(n_angles_);
    for (int pos = 0; pos < n_angles_; ++pos) degree_rank_[degree_order_[pos]] = pos;

    degree_to_idx_.assign(360, -1);
    for (int i = 0; i < n_angles_; ++i)
        degree_to_idx_[((bao_.angleDegree(i) % 360) + 360) % 360] = i;
}

// Mueve n_moves ángulos activos DISTINTOS de nb, cada uno una magnitud
// aleatoria en (base_step, 2*base_step), reintentando ante colisión con otro
// ángulo activo. Factorizado de AngleShiftMultiPerturbation::perturb() para
// que AdaptiveAngleShiftPerturbation::perturb() reutilice exactamente la
// misma mecánica de movimiento con n_moves/base_step variables en vez de
// fijos. No evalúa la solución (responsabilidad del llamador).
static void applyAngleShiftMultiMove(BaoSolution* nb, int n_angles,
                                      int n_moves,
                                      const std::vector<int>& degree_order,
                                      const std::vector<int>& degree_rank,
                                      int base_step,
                                      const std::vector<int>* magnitudes_deg = nullptr,
                                      const BaoProblem* bao = nullptr,
                                      const std::vector<int>* degree_to_idx = nullptr,
                                      const std::function<void()>* after_move = nullptr)
{
    int K = (int)nb->active_angles_.size();
    n_moves = std::min(n_moves, K);   // no se puede mover más ángulos distintos que K

    // Permutación de slots (Fisher-Yates): garantiza n_moves slots DISTINTOS,
    // a diferencia de n_moves sorteos independientes con reposición que
    // pueden repetir el mismo slot (y terminar moviendo un solo ángulo en
    // vez de varios — el bug que se detectó revisando datos reales).
    std::vector<int> slot_order(K);
    for (int i = 0; i < K; ++i) slot_order[i] = i;
    for (int i = K - 1; i > 0; --i) {
        int j = emili::generateRandomNumber() % (i + 1);
        std::swap(slot_order[i], slot_order[j]);
    }

    // Magnitud aleatoria en (base_step, 2*base_step): igual razón que en
    // AngleShiftNeighborhood::random — nunca base_step exacto, para no
    // quedar atrapado en la misma clase módulo base_step que usa la
    // búsqueda local.
    base_step = ((base_step % n_angles) + n_angles) % n_angles;
    if (base_step < 1) base_step = 1;
    int span = base_step - 1;
    if (span < 1) span = 1;

    std::vector<bool> active_flag(n_angles, false);
    for (int a : nb->active_angles_) active_flag[a] = true;

    for (int m = 0; m < n_moves; ++m) {
        int ai = slot_order[m];
        const int max_attempts = 8;
        for (int attempt = 0; attempt < max_attempts; ++attempt) {
            if (magnitudes_deg) {
                // Modo magnitud discreta: +-m grados reales (no posiciones del
                // catalogo). Si el grado destino no existe o esta ocupado,
                // reintenta; si se agotan los intentos el slot queda igual.
                int mag = (*magnitudes_deg)[emili::generateRandomNumber() % magnitudes_deg->size()];
                int sign = (emili::generateRandomNumber() % 2) ? 1 : -1;
                int cur_idx = nb->active_angles_[ai];
                int tgt_deg = (((bao->angleDegree(cur_idx) + sign * mag) % 360) + 360) % 360;
                int cand = (*degree_to_idx)[tgt_deg];
                if (cand < 0 || active_flag[cand]) continue;
                active_flag[cur_idx] = false;
                active_flag[cand] = true;
                nb->active_angles_[ai] = cand;
                break;
            }
            int dir = emili::generateRandomNumber() % 2;
            int step_mod = (base_step > 1)
                ? base_step + 1 + (emili::generateRandomNumber() % span)
                : 1 + (emili::generateRandomNumber() % (n_angles - 1));

            int active_catalog_idx = nb->active_angles_[ai];
            int pos = degree_rank[active_catalog_idx];
            int new_pos = (dir == 0)
                ? (pos - step_mod + n_angles) % n_angles
                : (pos + step_mod) % n_angles;
            int candidate = degree_order[new_pos];

            if (active_flag[candidate]) continue;   // colisión: reintenta

            active_flag[active_catalog_idx] = false;
            active_flag[candidate] = true;
            nb->active_angles_[ai] = candidate;
            break;
        }
        // Modo stepwise: el llamador evalua el estado intermedio tras cada slot
        // (aunque no se haya podido mover). No consume numeros aleatorios.
        if (after_move) (*after_move)();
    }
}

emili::Solution* AngleShiftMultiPerturbation::perturb(emili::Solution* solution)
{
    if ((int)degree_order_.size() != n_angles_) buildDegreeOrder();

    BaoSolution* nb = static_cast<BaoSolution*>(solution->clone());
    int K = (int)nb->active_angles_.size();
    int n_moves = std::min(numSteps_, K);   // no se puede mover más ángulos distintos que K

    std::function<void()> eval_each = [this, nb]() { bao_.evaluateSolution(*nb); };
    applyAngleShiftMultiMove(nb, n_angles_, n_moves, degree_order_, degree_rank_, step_,
                              magnitudes_.empty() ? nullptr : &magnitudes_, &bao_, &degree_to_idx_,
                              stepwise_ ? &eval_each : nullptr);
    if (stepwise_) return nb;   // el ultimo movimiento ya evaluo el estado final

    // Por defecto: una sola resolución FMO al final, con los n_moves ángulos
    // ya desplazados — los estados intermedios nunca se usan para decidir
    // nada, evaluarlos aparte solo gasta resoluciones de más. Opt-in
    // (stepwise_, token "stepwise"): se evalúa tras CADA movimiento
    // individual (n_moves resoluciones, cada estado intermedio queda
    // registrado como cualquier otra evaluación); la solución devuelta sigue
    // siendo el estado final y la secuencia aleatoria es idéntica.
    bao_.evaluateSolution(*nb);
    return nb;
}

void AdaptiveAngleShiftPerturbation::buildDegreeOrder()
{
    // Idéntico a AngleShiftMultiPerturbation::buildDegreeOrder — ver ese comentario.
    degree_order_.resize(n_angles_);
    for (int i = 0; i < n_angles_; ++i) degree_order_[i] = i;
    std::sort(degree_order_.begin(), degree_order_.end(),
        [this](int a, int b) { return bao_.angleDegree(a) < bao_.angleDegree(b); });

    degree_rank_.resize(n_angles_);
    for (int pos = 0; pos < n_angles_; ++pos) degree_rank_[degree_order_[pos]] = pos;

    degree_to_idx_.assign(360, -1);
    for (int i = 0; i < n_angles_; ++i)
        degree_to_idx_[((bao_.angleDegree(i) % 360) + 360) % 360] = i;
}

emili::Solution* AdaptiveAngleShiftPerturbation::perturb(emili::Solution* solution)
{
    if ((int)degree_order_.size() != n_angles_) buildDegreeOrder();

    // Auto-observación: compara el objetivo recibido en esta llamada contra
    // el recibido en la llamada anterior (perturb() ya recibe la solución
    // actual aceptada por el ILS en cada ronda -- ver
    // emili::IteratedLocalSearch::search, "s_p = pert.perturb(s);"). No hace
    // falta ningún hook nuevo en el framework genérico.
    double cur_obj = solution->getSolutionValue();
    if (cur_obj < last_objective_) {
        // Mejora real desde la última llamada: vuelve a nivel 0 (1 ángulo).
        level_ = 0;
        stagnation_count_ = 0;
    } else {
        stagnation_count_++;
        if (stagnation_count_ >= stagnationThreshold_) {
            // Nivel máximo: el que hace level_+1 == maxNumSteps_.
            int max_level = std::max(0, maxNumSteps_ - 1);
            level_ = std::min(level_ + 1, max_level);
            stagnation_count_ = 0;
        }
    }
    last_objective_ = cur_obj;

    // Nivel 0 mueve 1 ángulo; cada nivel adicional suma uno más, hasta
    // maxNumSteps_. La magnitud (step_) se mantiene fija en todos los
    // niveles -- solo escala la cantidad de ángulos movidos.
    int effective_numSteps = std::min(level_ + 1, maxNumSteps_);

    BaoSolution* nb = static_cast<BaoSolution*>(solution->clone());
    std::function<void()> eval_each = [this, nb]() { bao_.evaluateSolution(*nb); };
    applyAngleShiftMultiMove(nb, n_angles_, effective_numSteps, degree_order_,
                              degree_rank_, step_,
                              magnitudes_.empty() ? nullptr : &magnitudes_, &bao_, &degree_to_idx_,
                              stepwise_ ? &eval_each : nullptr);

    if (!stepwise_) bao_.evaluateSolution(*nb);
    return nb;
}

//=== Perturbación =============================================================

emili::Solution* RandomAnglesPerturbation::perturb(emili::Solution* current)
{
    // Trabaja sobre una copia — no muta la solución actual
    BaoSolution* bs = static_cast<BaoSolution*>(current->clone());
    const int n = bao_.nAngles();

    // Construye la lista de ángulos inactivos
    std::vector<bool> active_flag(n, false);
    for (int ai : bs->active_angles_) active_flag[ai] = true;
    std::vector<int> inactive;
    inactive.reserve(n - (int)bs->active_angles_.size());
    for (int a = 0; a < n; ++a)
        if (!active_flag[a]) inactive.push_back(a);

    // Limita los swaps al mínimo entre p_, activos disponibles e inactivos disponibles
    int swaps = std::min(p_, std::min((int)bs->active_angles_.size(),
                                      (int)inactive.size()));

    // Fisher-Yates parcial sobre ambas listas para seleccionar 'swaps' pares aleatorios
    for (int i = 0; i < swaps; ++i) {
        int ai = i + (int)(emili::generateRandomNumber() % (bs->active_angles_.size() - i));
        std::swap(bs->active_angles_[i], bs->active_angles_[ai]);
        int ii = i + (int)(emili::generateRandomNumber() % (inactive.size() - i));
        std::swap(inactive[i], inactive[ii]);
    }

    // Reemplaza los primeros 'swaps' activos por los primeros 'swaps' inactivos seleccionados
    for (int i = 0; i < swaps; ++i)
        bs->active_angles_[i] = inactive[i];

    std::sort(bs->active_angles_.begin(), bs->active_angles_.end());
    bao_.evaluateSolution(*bs);
    return bs;
}

//=== Aceptación ================================================================

emili::Solution* BaoImproveAccept::accept(emili::Solution* current,
                                           emili::Solution* candidate)
{
    bool improves = (*candidate < *current);
    if (!reject_repeated_) {
        return improves ? candidate : current;
    }

    // current siempre queda registrado, sin importar si mejora — así una
    // solución ya visitada se reconoce aunque se vuelva a alcanzar como
    // "current" en una iteración posterior.
    const BaoSolution* bc_cur = static_cast<const BaoSolution*>(current);
    visited_.insert(bc_cur->active_angles_);

    const BaoSolution* bc_cand = static_cast<const BaoSolution*>(candidate);
    if (improves && visited_.count(bc_cand->active_angles_)) {
        // candidate repite un conjunto de ángulos ya explorado — se rechaza
        // aunque mejore, para que ILS siga diversificando en vez de volver
        // a asentarse en un óptimo local ya conocido. Se imprime siempre
        // (no solo en verbose_) porque ocurre a lo sumo una vez por
        // iteración externa del ILS, no por cada evaluación FMO.
        std::cout << "  BAO accept: candidate f=" << std::fixed << std::setprecision(2)
                   << candidate->getSolutionValue()
                   << " repeats a previously visited angle set -- rejected, keeping current\n";
        return current;
    }

    visited_.insert(bc_cand->active_angles_);
    return improves ? candidate : current;
}

emili::Solution* MultiScaleAngleShake::shake(emili::Solution* s, int k)
{
    // k controla la intensidad del shake: k=0 → 1 swap, k=1 → 2 swaps, etc.
    // Delega en RandomAnglesPerturbation con p = k+1
    RandomAnglesPerturbation pert(bao_, k + 1);
    return pert.perturb(s);
}

emili::Solution* GreedyAnglesPerturbation::perturb(emili::Solution* current)
{
    // Trabaja sobre una copia — no muta la solución actual
    BaoSolution* bs = static_cast<BaoSolution*>(current->clone());
    const int n = bao_.nAngles();
    const int K = bao_.K();

    // Construye la lista de ángulos inactivos
    std::vector<bool> active_flag(n, false);
    for (int ai : bs->active_angles_) active_flag[ai] = true;
    std::vector<int> inactive;
    inactive.reserve(n - K);
    for (int a = 0; a < n; ++a)
        if (!active_flag[a]) inactive.push_back(a);

    // D: cuántos ángulos destruir y reconstruir — acotado por activos e inactivos disponibles
    int D = std::min(D_, std::min(K, (int)inactive.size()));

    // Destrucción: elimina D ángulos activos al azar (Fisher-Yates sobre el prefijo)
    for (int i = 0; i < D; ++i) {
        int ai = i + emili::generateRandomNumber() % (K - i);
        std::swap(bs->active_angles_[i], bs->active_angles_[ai]);
        inactive.push_back(bs->active_angles_[i]);   // el ángulo removido pasa a inactivos
    }
    bs->active_angles_.erase(bs->active_angles_.begin(),
                              bs->active_angles_.begin() + D);

    // Construcción greedy: agrega D ángulos uno a uno eligiendo el mejor FMO en cada paso
    for (int step = 0; step < D; ++step) {
        double best_f   = 1e30;
        int    best_pos = 0;

        for (int j = 0; j < (int)inactive.size(); ++j) {
            // Prueba agregar inactive[j] temporalmente y evalúa
            bs->active_angles_.push_back(inactive[j]);
            std::sort(bs->active_angles_.begin(), bs->active_angles_.end());

            double f = bao_.calcObjectiveFunctionValue(*bs);
            bs->setSolutionValue(f);

            if (f < best_f) { best_f = f; best_pos = j; }

            // Deshace: elimina el ángulo de prueba
            bs->active_angles_.erase(
                std::find(bs->active_angles_.begin(),
                          bs->active_angles_.end(), inactive[j]));
        }

        // Confirma el mejor ángulo encontrado en este paso
        bs->active_angles_.push_back(inactive[best_pos]);
        std::sort(bs->active_angles_.begin(), bs->active_angles_.end());
        bao_.evaluateSolution(*bs);
        inactive.erase(inactive.begin() + best_pos);
    }

    return bs;
}

//=== Memoria tabú =============================================================

bool BaoTabuMemory::tabu_check(emili::Solution* s)
{
    const BaoSolution* bs = static_cast<const BaoSolution*>(s);

    // Recorre el buffer circular desde la entrada más reciente hacia atrás
    for (int i = 0; i < count_; ++i) {
        int idx = (head_ - 1 - i + tenure_) % tenure_;   // índice circular
        if (memory_[idx] == bs->active_angles_)
            return true;   // la solución ya fue visitada — está en la lista tabú
    }
    return false;
}

void BaoTabuMemory::forbid(emili::Solution* s)
{
    const BaoSolution* bs = static_cast<const BaoSolution*>(s);

    // Escribe en la posición actual del buffer y avanza el puntero circular
    memory_[head_] = bs->active_angles_;
    head_ = (head_ + 1) % tenure_;

    // count_ crece hasta tenure_ y se mantiene ahí — el buffer es de tamaño fijo
    if (count_ < tenure_) ++count_;
}

void AdaptiveBaoTabuMemory::resizeTenure(int new_tenure)
{
    // Crea un nuevo buffer con el nuevo tamaño
    std::vector<std::vector<int>> new_mem(new_tenure);

    // Copia las entradas más recientes que quepan en el nuevo buffer
    int entries = std::min(count_, new_tenure);
    for (int i = 0; i < entries; ++i) {
        int old_idx = (head_ - 1 - i + tenure_) % tenure_;
        int new_idx = (new_tenure - 1 - i + new_tenure) % new_tenure;
        new_mem[new_idx] = memory_[old_idx];
    }

    // Reemplaza el buffer y actualiza el estado interno
    memory_  = std::move(new_mem);
    head_    = entries % new_tenure;
    count_   = entries;
    tenure_  = new_tenure;
    setTabuTenure(new_tenure);
}

void AdaptiveBaoTabuMemory::forbid(emili::Solution* s)
{
    if (tabu_check(s)) {
        // Oscilación detectada: la solución ya fue visitada → aumenta el tenure para escapar
        if (tenure_ < tenure_max_) resizeTenure(tenure_ + 1);
        since_last_revisit_ = 0;
    } else {
        ++since_last_revisit_;
        // Período tranquilo: sin revisitas por 2×tenure pasos → reduce tenure para explorar más
        if (since_last_revisit_ >= 2 * tenure_ && tenure_ > tenure_min_) {
            resizeTenure(tenure_ - 1);
            since_last_revisit_ = 0;
        }
    }

    // Registra la solución en el buffer tabú base
    BaoTabuMemory::forbid(s);
}

} // namespace imrt
} // namespace emili
