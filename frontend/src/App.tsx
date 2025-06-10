import Layout from "./components/Layout";
import LoginForm from "./components/LoginForm";
import RegisterForm from "./components/RegisterForm";
import VerifyPage from "./components/VerifyPage";
import ProtectedRoute from "./components/ProtectedRoute";
import CreateEnvironment from "./components/CreateEnvironment";
import Dashboard from "./components/Dashboard";
import SettingsEnvironments from "./components/SettingsEnvironments";
import TermsPrivacy from "./components/TermsPrivacy";
import Settings from "./components/Settings";
import { BrowserRouter as Router, Routes, Route } from "react-router-dom";

function App() {
  return (
    <div className="w-screen h-screen min-h-screen min-w-full">
      <Router>
        <Routes>
          {/* Public routes: no Layout */}
          <Route path="/login" element={<LoginForm />} />
          <Route path="/register" element={<RegisterForm />} />
          <Route path="/verify" element={<VerifyPage />} />
          <Route path="/policies/terms-of-use/" element={<TermsPrivacy />} />
          <Route path="/policies/privacy-policy/" element={<TermsPrivacy />} />

          {/* Protected routes: with Layout */}
          <Route
            path="/"
            element={
              <ProtectedRoute>
                <Layout>
                  <Dashboard />
                </Layout>
              </ProtectedRoute>
            }
          />
          <Route
            path="/environments/new"
            element={
              <ProtectedRoute>
                <Layout>
                  <CreateEnvironment />
                </Layout>
              </ProtectedRoute>
            }
          />
          <Route
            path="/settings"
            element={
              <ProtectedRoute>
                <Layout>
                  <Settings />
                </Layout>
              </ProtectedRoute>
            }
          />
          <Route
            path="/settings/environments"
            element={
              <ProtectedRoute>
                <Layout>
                  <SettingsEnvironments />
                </Layout>
              </ProtectedRoute>
            }
          />
        </Routes>
      </Router>
    </div>
  );
}

export default App;
