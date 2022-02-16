#pragma once

#include <Hull/Hull.h>
#include <math.h>
#include <memory>
#include <string>

namespace hull {
void toObj(const std::vector<Coordinate> &vertices,
           const std::vector<Facet> &facets, const std::string &fileName);

bool check_normals(const std::vector<Coordinate> &vertices,
                   const std::vector<Facet> &facets);

class TrivialObserver : public Observer {
public:
  TrivialObserver() = default;

  std::unique_ptr<Notification> last_notification;
  void hullChanges(const Notification &notification) override {
    last_notification = std::make_unique<Notification>(notification);
  };
};

class StepsLogger : public hull::Observer {
public:
  StepsLogger(const std::string &log_name) : log_name(log_name){};

  void hullChanges(const Notification &notification) override;

  void check_updated_mesh(const Notification &notification) const;

private:
  const std::string log_name;
};
} // namespace hull

std::string generate_obj_log_name(const std::string &file_name);

constexpr float GREEK_PI = 3.14159f;
float to_rad(const float angle);

std::vector<float> linspace(const float min, const float max,
                            const std::size_t intervals);

hull::Hull fill_hull(const std::vector<hull::Coordinate> &vertices,
                     hull::Observer *obs = nullptr);

bool are_same(const std::vector<hull::Coordinate> &a,
              const std::vector<hull::Coordinate> &b);
