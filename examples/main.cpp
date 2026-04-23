#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "../c_ffi.hpp"

int main()
{
    std::vector<vec2> polygon = {
        {0.0, 0.0},
        {0.0, 30.0},
        {30.0, 30.0},
        {30.0, 0.0},
        {0.0, 0.0},
    };
    
    PhaseGrid *const phase_grid = PhaseGrid_new(
        polygon.data(),
        polygon.size(),
        {10.0, 10.0}
    );

    std::vector<LiquidDropState> initial_states = {
        {
            {0.5, 1.0},
            {1.0, 1.0},
            std::pow(0.1, 3.0)
        },
    };

    Rk4 *rk4 = Rk4_new(
        1.2041,
        1.8e-05,
        1000.0,
        72.8,
        500.0,
        10.0,
        0.99,
        1.0,
        phase_grid,
        initial_states.data(),
        initial_states.size(),
        0.0,
        1.0,
        0.1
    );
    
    Rk4Result *result = Rk4_integrate(rk4);
    auto sequences_len = Rk4Result_len(result);
    auto sequences_data = Rk4Result_data(result);

    std::cout << "Общее число капель: " << sequences_len << std::endl;
    for (std::size_t i = 0; i < sequences_len; ++i) {
        std::cout << "Капля #" << i + 1 << std::endl;
        std::cout << "Число состояний: " << sequences_data[i].states_len << std::endl;
        std::cout << "Начало: " << sequences_data[i].time_begin << std::endl;
        std::cout << "Конец: " << sequences_data[i].time_end << std::endl;
        std::cout << "Конец? " << (sequences_data[i].is_alive ? "нет" : "да") << std::endl;
        for (std::size_t j = 0; j < sequences_data[i].states_len; ++j) {
            std::cout << "position={"
                          << sequences_data[i].states_data[j].position.x << ","
                          << sequences_data[i].states_data[j].position.y << "}, "
                      << "speed={"
                          << sequences_data[i].states_data[j].speed.x << ","
                          << sequences_data[i].states_data[j].speed.y << "}, "
                      << "diameter=" << sequences_data[i].states_data[j].diameter3 << ", "
                      << "accumulated_stress=" << sequences_data[i].states_data[j].accumulated_stress
                      << std::endl;
            
        }
    }

    Rk4Result_destroy(result);
    Rk4_destroy(rk4);
    PhaseGrid_destroy(phase_grid);
    return 0;
}