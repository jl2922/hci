#ifndef CONNECTIONS_H_
#define CONNECTIONS_H_

class Connections {
 public:
  virtual ~Connections() = default;

  virtual void update() = 0;

  virtual std::vector<std::pair<int, double>> get_conections(const int i) = 0;
};

#endif
