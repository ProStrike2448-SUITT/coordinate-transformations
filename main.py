import inspect
import math
import random
import time

import matplotlib.pyplot as plt

# Типи даних
PolarPoint = tuple[float, float]  # (r, theta)
CartesianPoint2D = tuple[float, float]  # (x, y)
SphericalPoint = tuple[float, float, float]  # (r, theta, phi)
CartesianPoint3D = tuple[float, float, float]  # (x, y, z)

# ------------------------------
# Функції для перетворення координат
# ------------------------------


def polar_to_cartesian(r: float, theta: float) -> CartesianPoint2D:
    """Перетворення з полярної системи в декартову."""
    x = r * math.cos(theta)
    y = r * math.sin(theta)
    return x, y


def cartesian_to_polar(x: float, y: float) -> PolarPoint:
    """Перетворення з декартової системи в полярну."""
    r = math.sqrt(x**2 + y**2)
    theta = math.atan2(y, x) % (2 * math.pi)  # Нормалізація в межах [0, 2π)
    return r, theta


def spherical_to_cartesian(r: float, theta: float, phi: float) -> CartesianPoint3D:
    """Перетворення з сферичної системи в декартову."""
    x = r * math.sin(phi) * math.cos(theta)
    y = r * math.sin(phi) * math.sin(theta)
    z = r * math.cos(phi)
    return x, y, z


def cartesian_to_spherical(x: float, y: float, z: float) -> SphericalPoint:
    """Перетворення з декартової системи в сферичну."""
    r = math.sqrt(x**2 + y**2 + z**2)
    theta = math.atan2(y, x) % (2 * math.pi)  # Нормалізація в межах [0, 2π)
    phi = math.acos(z / r) if r != 0 else 0
    return r, theta, phi


# ------------------------------
# Функції для обчислення відстаней
# ------------------------------


def distance_cartesian_2d(x1: float, y1: float, x2: float, y2: float) -> float:
    """Пряма відстань у 2D в декартовій системі."""
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


def distance_polar(r1: float, theta1: float, r2: float, theta2: float) -> float:
    """Відстань у полярній системі."""
    return math.sqrt(r1**2 + r2**2 - 2 * r1 * r2 * math.cos(theta2 - theta1))


def distance_cartesian_3d(
    x1: float, y1: float, z1: float, x2: float, y2: float, z2: float
) -> float:
    """Пряма відстань у 3D в декартовій системі."""
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)


def distance_spherical_volume(
    r1: float, theta1: float, phi1: float, r2: float, theta2: float, phi2: float
) -> float:
    """Відстань через об'єм сфери."""
    x1, y1, z1 = spherical_to_cartesian(r1, theta1, phi1)
    x2, y2, z2 = spherical_to_cartesian(r2, theta2, phi2)
    return distance_cartesian_3d(x1, y1, z1, x2, y2, z2)


def distance_spherical_surface(
    r: float, theta1: float, phi1: float, theta2: float, phi2: float
) -> float:
    """Дугова відстань по поверхні сфери."""
    return r * math.acos(
        math.sin(phi1) * math.sin(phi2)
        + math.cos(phi1) * math.cos(phi2) * math.cos(theta2 - theta1)
    )


# ------------------------------
# Генерація випадкових координат
# ------------------------------


def generate_random_polar_points(n: int) -> list[PolarPoint]:
    """Генерує випадкові точки у полярній системі."""
    return [(random.uniform(0, 100), random.uniform(0, 2 * math.pi)) for _ in range(n)]


def generate_random_spherical_points(n: int) -> list[SphericalPoint]:
    """Генерує випадкові точки у сферичній системі."""
    return [
        (random.uniform(0, 100), random.uniform(0, 2 * math.pi), random.uniform(0, math.pi))
        for _ in range(n)
    ]


# ------------------------------
# Бенчмаркінг
# ------------------------------


def benchmark(func, points: list[tuple[float, ...]]) -> float:
    """Вимірює час виконання функції для масиву точок."""
    start_time = time.perf_counter()
    if len(inspect.signature(func).parameters) == 5:
        [func(*p1, *p2[1:]) for p1, p2 in zip(points[::2], points[1::2])]
    else:
        [func(*p1, *p2) for p1, p2 in zip(points[::2], points[1::2])]
    end_time = time.perf_counter()
    return end_time - start_time


# ------------------------------
# Головна програма
# ------------------------------

if __name__ == "__main__":
    random.seed(214215)
    num_points = 100000

    # Двовимірний простір: Полярна -> Декартова -> Полярна
    polar_points = generate_random_polar_points(num_points)
    cartesian_points = [polar_to_cartesian(r, theta) for r, theta in polar_points]
    converted_polar_points = [cartesian_to_polar(x, y) for x, y in cartesian_points]

    # Перевірка коректності
    assert all(
        math.isclose(p1[0], p2[0], rel_tol=1e-5) and math.isclose(p1[1], p2[1], rel_tol=1e-5)
        for p1, p2 in zip(polar_points, converted_polar_points)
    ), "Перетворення некоректне!"

    # Тривимірний простір: Сферична -> Декартова -> Сферична
    spherical_points = generate_random_spherical_points(num_points)
    cartesian_points_3d = [
        spherical_to_cartesian(r, theta, phi) for r, theta, phi in spherical_points
    ]
    converted_spherical_points = [
        cartesian_to_spherical(x, y, z) for x, y, z in cartesian_points_3d
    ]

    # Перевірка коректності
    assert all(
        math.isclose(p1[0], p2[0], rel_tol=1e-5)
        and math.isclose(p1[1], p2[1], rel_tol=1e-5)
        and math.isclose(p1[2], p2[2], rel_tol=1e-5)
        for p1, p2 in zip(spherical_points, converted_spherical_points)
    ), "Перетворення некоректне!"

    # Бенчмаркінг
    polar_time = benchmark(distance_polar, polar_points)
    cartesian_time_2d = benchmark(distance_cartesian_2d, cartesian_points)
    cartesian_time_3d = benchmark(distance_cartesian_3d, cartesian_points_3d)
    spherical_time_volume = benchmark(
        distance_spherical_volume,
        spherical_points,
    )
    spherical_time_surface = benchmark(
        distance_spherical_surface,
        spherical_points,
    )

    coordinate_systems = [
        "Полярна",
        "Декартова\n2D",
        "Декартова\n3D",
        "Сферична\n(об'єм)",
        "Сферична\n(поверхня)",
    ]
    times = [
        polar_time,
        cartesian_time_2d,
        cartesian_time_3d,
        spherical_time_volume,
        spherical_time_surface,
    ]

    # Графік
    fig, ax = plt.subplots(figsize=(6, 7))
    p = ax.bar(
        coordinate_systems,
        times,
        color=["blue", "green", "red", "orange", "purple"],
    )
    ax.bar_label(p, [round(time, 6) for time in times])

    ax.set_title("Порівняння часу необхідного для обчислення відстані\nу різних системах координат")
    ax.set_ylabel("Час виконання (с)")
    ax.set_xlabel("Система координат")

    plt.show()
