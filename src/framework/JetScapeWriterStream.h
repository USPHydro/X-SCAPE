/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 *
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

// jetscape writer ascii class

#ifndef JETSCAPEWRITERSTREAM_H
#define JETSCAPEWRITERSTREAM_H

#include <fstream>
#include <string>

#ifdef USE_GZIP
#include "gzstream.h"
#endif

#include "JetScapeWriter.h"

using std::ofstream;

namespace Jetscape {

template <class T> class JetScapeWriterStream : public JetScapeWriter {

public:
  JetScapeWriterStream<T>(){};
  JetScapeWriterStream<T>(string m_file_name_out);
  virtual ~JetScapeWriterStream<T>();

  void InitTask();
  void ExecuteTask();

  bool GetStatus() { return output_file.good(); }
  void Close() { output_file.close(); }

  void WriteInitFileXMLMain();
  void WriteInitFileXMLUser();

  void Write(weak_ptr<PartonShower> ps);
  void Write(weak_ptr<Parton> p);
  void Write(weak_ptr<Vertex> v);
  void Write(weak_ptr<Hadron> h);
  void WriteHeaderToFile();

  void Write(string s) { output_file << s << endl; }
  void WriteComment(string s) { output_file << "# " << s << endl; }
  void WriteWhiteSpace(string s) { output_file << s << " "; }
  void WriteEvent();

protected:
  T output_file; //!< Output file
  int hadronCounter = 0;
  //int m_precision; //!< Output precision

  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<JetScapeWriterStream<ofstream>> reg;
  static RegisterJetScapeModule<JetScapeWriterStream<ogzstream>> regGZ;
};

typedef JetScapeWriterStream<ofstream> JetScapeWriterAscii;
#ifdef USE_GZIP
typedef JetScapeWriterStream<ogzstream> JetScapeWriterAsciiGZ;
#endif

} // end namespace Jetscape

#endif // JETSCAPEWRITERSTREAM_H
