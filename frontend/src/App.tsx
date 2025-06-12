import Layout from "./components/Layout";
import LoginForm from "./components/LoginForm";
import RegisterForm from "./components/RegisterForm";
import VerifyPage from "./components/VerifyPage";
import ProtectedRoute from "./components/ProtectedRoute";
import CreateEnvironment from "./components/CreateEnvironment";
import Dashboard from "./components/Dashboard";
import SettingsEnvironments from "./components/SettingsEnvironments";
import TermsPrivacy from "./components/TermsPrivacy";
import { SettingsGeneral } from "./components/Settings";
import SettingsDataControls from "./components/SettingsDataControls";
import CreateEnvironmentSettings from "./components/CreateEnvironmentSettings";
import { BrowserRouter as Router, Routes, Route } from "react-router-dom";
import { EnvironmentDetails } from "./components/CreateEnvironment";
import { useState } from "react";
import ChatSession from "./components/ChatSession";

function App() {
  const [selectedModel, setSelectedModel] = useState<string>("vanilla");
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
                  <Dashboard selectedModel={selectedModel} setSelectedModel={setSelectedModel} />
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
            path="/settings/general"
            element={
              <ProtectedRoute>
                <Layout>
                  <SettingsGeneral />
                </Layout>
              </ProtectedRoute>
            }
          />
          <Route
            path="/settings/data-controls"
            element={
              <ProtectedRoute>
                <Layout>
                  <SettingsDataControls />
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
          <Route
            path="/settings/environments/create"
            element={
              <ProtectedRoute>
                <Layout>
                  <CreateEnvironmentSettings />
                </Layout>
              </ProtectedRoute>
            }
          />
          <Route
            path="/environment/:env_id"
            element={
              <ProtectedRoute>
                <Layout>
                  <EnvironmentDetails />
                </Layout>
              </ProtectedRoute>
            }
          />
          <Route
            path="/chat/session/:sessionId"
            element={
              <ProtectedRoute>
                <Layout>
                  <ChatSession selectedModel={selectedModel} setSelectedModel={setSelectedModel} />
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
