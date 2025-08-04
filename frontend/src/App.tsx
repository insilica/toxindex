import Layout from "./components/Layout";
import LoginForm from "./components/LoginForm";
import RegisterForm from "./components/RegisterForm";
import VerifyPage from "./components/VerifyPage";
import ProtectedRoute from "./components/ProtectedRoute";
// import CreateEnvironment from "./components/CreateEnvironment";
import Dashboard from "./components/Dashboard";
import SettingsEnvironments from "./components/SettingsEnvironments";
import TermsPrivacy from "./components/TermsPrivacy";
import Settings from "./components/Settings";
import SettingsDataControls from "./components/SettingsDataControls";
import CreateEnvironmentSettings from "./components/CreateEnvironmentSettings";
import { BrowserRouter as Router, Routes, Route, useNavigate } from "react-router-dom";
import EnvironmentDetails from "./components/EnvironmentDetails";

import ChatSession from "./components/ChatSession";
import TaskDetail from './components/TaskDetail';
import UserProfile from './components/UserProfile';
import ForgotPasswordPage from './components/ForgotPasswordPage';
import ResetPasswordPage from './components/ResetPasswordPage';
import AdminRoute from './components/admin/AdminRoute';
import UserGroupManager from './components/admin/UserGroupManager';
import { EnvironmentProvider } from "./context/EnvironmentContext";
import { ChatSessionProvider } from "./context/ChatSessionContext";
import { ModelProvider } from "./context/ModelContext";
import { SocketProvider } from "./context/SocketContext";
import { SessionProvider, useSession } from "./context/SessionContext";
import { SessionTimeoutNotification } from "./components/SessionTimeoutNotification";

function AppContent() {
  const navigate = useNavigate();

  const handleSessionExpired = () => {
    console.log('[App] Session expired, redirecting to login');
    // Clear any stored auth data
    localStorage.removeItem('authToken');
    sessionStorage.clear();
    // Redirect to login
    navigate('/login');
  };

  // Get session state from context
  const {
    extendSession,
    logoutNow,
    showWarning,
    timeRemaining
  } = useSession();

  return (
    <>
      <SessionTimeoutNotification
        isVisible={showWarning}
        onExtendSession={extendSession}
        onLogout={logoutNow}
        timeRemaining={timeRemaining}
      />
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
                <Dashboard />
              </Layout>
            </ProtectedRoute>
          }
        />
        <Route
          path="/settings/general"
          element={
            <ProtectedRoute>
              <Layout>
                <Settings/>
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
          path="/environments/details"
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
                <ChatSession />
              </Layout>
            </ProtectedRoute>
          }
        />
        <Route path="/task/:task_id" element={<TaskDetail />} />
        <Route path="/user/:user_id" element={<UserProfile />} />

        {/* Admin routes */}
        <Route
          path="/admin/users"
          element={
            <ProtectedRoute>
              <Layout>
                <AdminRoute>
                  <UserGroupManager />
                </AdminRoute>
              </Layout>
            </ProtectedRoute>
          }
        />
      </Routes>
    </>
  );
}

function App() {
  console.log("App mounted");
  
  const handleSessionExpired = () => {
    // Clear any cached data
    localStorage.clear();
    sessionStorage.clear();
    // The navigation will be handled by AppContent
    window.location.href = '/login';
  };

  return (
    <div className="w-screen h-screen min-h-screen min-w-full">
      <SocketProvider>
        <ModelProvider>
          <ChatSessionProvider>
            <EnvironmentProvider>
              <Router>
                <SessionProvider onSessionExpired={handleSessionExpired}>
                  <AppContent />
                </SessionProvider>
              </Router>
            </EnvironmentProvider>
          </ChatSessionProvider>
        </ModelProvider>
      </SocketProvider>
    </div>
  );
}

export default App;
