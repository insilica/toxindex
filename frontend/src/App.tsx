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
import EnvironmentDetails from "./components/EnvironmentDetails";
import { useState, useCallback } from "react";
import ChatSession from "./components/ChatSession";
import TaskDetail from './components/TaskDetail';
import UserProfile from './components/UserProfile';
import ForgotPasswordPage from './components/ForgotPasswordPage';
import ResetPasswordPage from './components/ResetPasswordPage';

function App() {
  const [environments, setEnvironments] = useState<any[]>([]);
  const [loadingEnvironments, setLoadingEnvironments] = useState<boolean>(false);

  // Robust refetch function for environments
  const refetchEnvironments = useCallback(() => {
    setLoadingEnvironments(true);
    return fetch('/api/environments', { credentials: 'include' })
      .then(res => res.json())
      .then(data => {
        setEnvironments(data.environments || []);
        return data.environments || [];
      })
      .finally(() => setLoadingEnvironments(false));
  }, []);

  // Function to refresh files for a specific environment
  const refreshEnvFiles = useCallback((envId: string) => {
    return fetch(`/api/environments/${envId}/files`, { 
      credentials: 'include',
      cache: 'no-store'
    })
      .then(res => res.json())
      .then(data => data.files || []);
  }, []);

  return (
    <div className="w-screen h-screen min-h-screen min-w-full">
      <Router>
        <Routes>
          {/* Public routes: no Layout */}
          <Route path="/login" element={<LoginForm />} />
          <Route path="/register" element={<RegisterForm />} />
          <Route path="/verify/:token" element={<VerifyPage />} />
          <Route path="/policies/terms-of-use/" element={<TermsPrivacy />} />
          <Route path="/policies/privacy-policy/" element={<TermsPrivacy />} />
          <Route path="/forgot_password" element={<ForgotPasswordPage />} />
          <Route path="/reset_password/:token" element={<ResetPasswordPage />} />

          {/* Protected routes: with Layout */}
          <Route
            path="/"
            element={
              <ProtectedRoute>
                <Layout>
                  <Dashboard
                    environments={environments}
                    refetchEnvironments={refetchEnvironments}
                    loadingEnvironments={loadingEnvironments}
                    refreshEnvFiles={refreshEnvFiles}
                  />
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
                  <SettingsEnvironments
                    environments={environments}
                    refetchEnvironments={refetchEnvironments}
                    setEnvironments={setEnvironments}
                  />
                </Layout>
              </ProtectedRoute>
            }
          />
          <Route
            path="/settings/environments/create"
            element={
              <ProtectedRoute>
                <Layout>
                  <CreateEnvironmentSettings
                    refetchEnvironments={refetchEnvironments}
                  />
                </Layout>
              </ProtectedRoute>
            }
          />
          <Route
            path="/environment/:env_id"
            element={
              <ProtectedRoute>
                <Layout>
                  <EnvironmentDetails 
                    environments={environments}
                    loadingEnvironments={loadingEnvironments}
                    refetchEnvironments={refetchEnvironments}
                  />
                </Layout>
              </ProtectedRoute>
            }
          />
          <Route
            path="/chat/session/:sessionId"
            element={
              <ProtectedRoute>
                <Layout>
                  <ChatSession
                    environments={environments}
                    loadingEnvironments={loadingEnvironments}
                  />
                </Layout>
              </ProtectedRoute>
            }
          />
          <Route path="/task/:task_id" element={<TaskDetail />} />
          <Route path="/user/:user_id" element={<UserProfile />} />
        </Routes>
      </Router>
    </div>
  );
}

export default App;
