#ifndef LOOPWEAVER_SETTINGS_H
#define LOOPWEAVER_SETTINGS_H

// Default settings read from a text file, to change them without
// recompiling. Each parameter is registered with its name and a reference to
// the variable holding its (compiled) default value. LoadOrCreate:
//  - if the file does not exist it is written with the current values;
//  - otherwise the values in the file replace the defaults, and the file is
//    written again only if some registered parameter is missing in it.
// File format, one parameter per line ('#' starts a comment):
//   name value   # description

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

class LoopWeaverSettings {

  enum EntryType { TBool, TInt, TFloat, TDouble };

  struct Entry {
    std::string Name;
    EntryType Type;
    void *Ptr;
    std::string Comment;
  };

  std::vector<Entry> Entries;

  void AddEntry(const std::string &Name, EntryType Type, void *Ptr,
                const std::string &Comment) {
    Entries.push_back({Name, Type, Ptr, Comment});
  }

  static std::string ValueString(const Entry &E) {
    std::ostringstream oss;
    switch (E.Type) {
    case TBool:
      oss << (*(bool *)E.Ptr ? 1 : 0);
      break;
    case TInt:
      oss << *(int *)E.Ptr;
      break;
    case TFloat:
      oss << *(float *)E.Ptr;
      break;
    case TDouble:
      oss << *(double *)E.Ptr;
      break;
    }
    return oss.str();
  }

  static bool ParseValue(const Entry &E, const std::string &Value) {
    std::istringstream iss(Value);
    switch (E.Type) {
    case TBool: {
      std::string V;
      iss >> V;
      if ((V == "1") || (V == "true"))
        *(bool *)E.Ptr = true;
      else if ((V == "0") || (V == "false"))
        *(bool *)E.Ptr = false;
      else
        return false;
      return true;
    }
    case TInt:
      return static_cast<bool>(iss >> *(int *)E.Ptr);
    case TFloat:
      return static_cast<bool>(iss >> *(float *)E.Ptr);
    case TDouble:
      return static_cast<bool>(iss >> *(double *)E.Ptr);
    }
    return false;
  }

public:
  // register a parameter (its current value is the default)
  void Add(const std::string &Name, bool &V, const std::string &Comment = "") {
    AddEntry(Name, TBool, &V, Comment);
  }
  void Add(const std::string &Name, int &V, const std::string &Comment = "") {
    AddEntry(Name, TInt, &V, Comment);
  }
  void Add(const std::string &Name, float &V, const std::string &Comment = "") {
    AddEntry(Name, TFloat, &V, Comment);
  }
  void Add(const std::string &Name, double &V,
           const std::string &Comment = "") {
    AddEntry(Name, TDouble, &V, Comment);
  }

  // a section title in the written file
  void Section(const std::string &Title) { AddEntry(Title, TBool, NULL, ""); }

  bool Save(const std::string &Path) const {
    std::ofstream out(Path);
    if (!out)
      return false;
    out << "# LoopWeaver default settings: name value  (# comment)\n";
    out << "# delete this file to write it again with the compiled defaults\n";
    for (size_t i = 0; i < Entries.size(); i++) {
      const Entry &E = Entries[i];
      if (E.Ptr == NULL) {
        out << "\n# --- " << E.Name << " ---\n";
        continue;
      }
      std::string Line = E.Name + " " + ValueString(E);
      out << Line;
      if (!E.Comment.empty()) {
        for (size_t k = Line.size(); k < 40; k++)
          out << " ";
        out << "  # " << E.Comment;
      }
      out << "\n";
    }
    return true;
  }

  // reads the file; returns false if it does not exist. NumMissing: the
  // registered parameters not in the file
  bool Load(const std::string &Path, int &NumMissing) {
    NumMissing = 0;
    std::ifstream in(Path);
    if (!in)
      return false;
    std::map<std::string, std::string> Values;
    std::string Line;
    int LineNum = 0;
    while (std::getline(in, Line)) {
      LineNum++;
      size_t Hash = Line.find('#');
      if (Hash != std::string::npos)
        Line = Line.substr(0, Hash);
      std::istringstream iss(Line);
      std::string Name, Value;
      if (!(iss >> Name))
        continue;
      if (!(iss >> Value)) {
        std::cout << "Settings: no value for " << Name << " (line " << LineNum
                  << ")" << std::endl;
        continue;
      }
      Values[Name] = Value;
    }
    for (size_t i = 0; i < Entries.size(); i++) {
      const Entry &E = Entries[i];
      if (E.Ptr == NULL)
        continue;
      auto It = Values.find(E.Name);
      if (It == Values.end()) {
        NumMissing++;
        continue;
      }
      if (!ParseValue(E, It->second))
        std::cout << "Settings: wrong value for " << E.Name << ": "
                  << It->second << " (default kept)" << std::endl;
      Values.erase(It);
    }
    for (auto &V : Values)
      std::cout << "Settings: unknown parameter " << V.first << " (ignored)"
                << std::endl;
    return true;
  }

  void LoadOrCreate(const std::string &Path) {
    int NumMissing = 0;
    if (!Load(Path, NumMissing)) {
      if (Save(Path))
        std::cout << "Settings: written the defaults on " << Path
                  << std::endl;
      else
        std::cout << "Settings: cannot write " << Path << std::endl;
      return;
    }
    std::cout << "Settings: loaded " << Path << std::endl;
    // add the new parameters to the file
    if (NumMissing > 0) {
      Save(Path);
      std::cout << "Settings: added " << NumMissing << " missing parameters to "
                << Path << std::endl;
    }
  }
};

#endif
