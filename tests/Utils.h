#pragma once

#include <Hull/Hull.h>

#include <algorithm>
#include <filesystem>
#include <math.h>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>

namespace hull {
void toObj(const HullContext &hull, const std::filesystem::path &fileName);

class Logger {
public:
  static Logger &get();

  std::filesystem::path makeLogfileName(const std::string &topicName);

private:
  Logger();

  std::unordered_map<std::string, std::size_t> topics;
};

bool check_normals(const HullContext &hull);

class TrivialObserver : public Observer {
public:
  TrivialObserver() = default;

  std::optional<Notification> last_notification;

  void hullChanges(const Notification &notification) override {
    last_notification.emplace(notification);
  };
};

class StepsLogger : public hull::Observer {
public:
  StepsLogger(const std::filesystem::path &log_name) : log_name(log_name){};

  void hullChanges(const Notification &notification) override;

  void check_updated_mesh(const Notification &notification) const;

private:
  const std::filesystem::path log_name;
};

} // namespace hull

static constexpr float GREEK_PI = 3.14159f;
float to_rad(const float angle);

std::vector<float> linspace(const float min, const float max,
                            std::size_t intervals);

hull::Hull fill_hull(const std::vector<hull::Coordinate> &vertices,
                     hull::Observer *obs = nullptr);

template <typename A, typename B> bool same_vertices(const A &a, const B &b) {
  if (a.size() != b.size()) {
    return false;
  }
  for (const hull::Coordinate &a_vertex : a) {
    if (std::find_if(b.begin(), b.end(),
                     [&a_vertex](const hull::Coordinate &b_vertex) {
                       hull::Coordinate diff;
                       hull::diff(diff, b_vertex, a_vertex);
                       return hull::norm(diff) < 1e-5;
                     }) == b.end()) {
    }
  }
  return true;
}